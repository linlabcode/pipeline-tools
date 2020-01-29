#!/usr/bin/python
"""
The MIT License (MIT)

Copyright (c) 2019 Charles Lin

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
"""

# ==========================================================================
# =============================DEPENDENCIES=================================
# ==========================================================================

import argparse
import os
from collections import defaultdict

import numpy
from pipeline_tools.utils import pipeline_utils, utils
from scipy import stats


def parse_args(args=None):
    """Argument parser."""
    parser = argparse.ArgumentParser(
        usage=(  # TODO: fix this section
            "out_degree_delta [options]"
            " --edge_table_1 [EDGE_TABLE_1]"
            " --edge_table_2 [EDGE_TABLE_2]"
            " --bams_1 [BAM1 BAM2 BAM3 ...]"
            " --bams_2 [BAM4 BAM5 BAM6 ...]"
            " -o [OUTPUTFOLDER]"
            " -n [NAME]"
        )
    )

    # Required flags
    parser.add_argument(
        "--edge_table_1",
        dest="edge_table_1",
        default=None,
        type=str,
        help=("Edge table for group 1"),
        required=True,
    )
    parser.add_argument(
        "--edge_table_2",
        dest="edge_table_2",
        default=None,
        type=str,
        help=("Edge table for group 2"),
        required=True,
    )
    parser.add_argument(
        "--bams_1",
        dest="bams_1",
        type=str,
        nargs="+",
        help=("List of bam files for group 1"),
        required=True,
    )
    parser.add_argument(
        "--bams_2",
        dest="bams_2",
        type=str,
        nargs="+",
        help=("List of bam files for group 2"),
        required=True,
    )

    parser.add_argument(
        "-o",
        "--output",
        dest="output",
        default=None,
        type=str,
        help="Enter an output folder",
        required=True,
    )
    parser.add_argument(
        "-n",
        "--name",
        dest="name",
        default=None,
        type=str,
        help="Provide a name for the project",
        required=True,
    )

    # Additional options
    parser.add_argument(
        "--names_1",
        dest="names_1",
        type=str,
        nargs="+",
        help=("List of names for group 1"),
        required=False,
    )
    parser.add_argument(
        "--names_2",
        dest="names_2",
        type=str,
        nargs="+",
        help=("List of names for group 2"),
        required=False,
    )

    return parser.parse_args(args)


def merge_edge_tables(edge_table_1, edge_table_2, output_table="merged_EDGE_TABLE.txt"):
    region_dict = defaultdict(list)
    new_regions = defaultdict(list)
    with open(edge_table_1, "r") as edge_table1, open(edge_table_2, "r") as edge_table2:
        header = edge_table1.readline()
        for line in edge_table1:
            source, target, chrs, Pstart, Pstop, name, score = line.split("\t")[:7]
            Pstart = int(Pstart)
            Pstop = int(Pstop)
            region_dict["_".join([source, target, chrs])].append(
                (Pstart, Pstop, name, score)
            )

        edge_table2.readline()
        for line in edge_table2:
            source, target, chrs, Pstart, Pstop, name, score = line.split("\t")[:7]
            Pstart = int(Pstart)
            Pstop = int(Pstop)
            region_dict["_".join([source, target, chrs])].append(
                (Pstart, Pstop, name, score)
            )

    for pos in region_dict:
        # Sorting will order the regions making them easier to work with
        region_dict[pos].sort()
        regions = region_dict[pos]
        while regions:
            i = 0
            region_1 = regions[i]
            start_1, stop_1, name_1, score_1 = region_1[0:4]
            while True:
                if len(regions) > i + 1:
                    region_2 = regions[i + 1]
                    start_2, stop_2 = region_2[0:2]
                else:
                    region_new = (str(start_1), str(stop_1), name_1, score_1)
                    break

                if start_2 <= stop_1:
                    stop_1 = stop_2
                    i += 1
                else:
                    region_new = (str(start_1), str(stop_1), name_1, score_1)
                    break

            new_regions[pos].append(region_new)
            del regions[: i + 1]

    with open(output_table, "w") as outfile:
        outfile.write(header)
        for pos in sorted(new_regions):
            for region in new_regions[pos]:
                outfile.write("\t".join(pos.split("_")) + "\t" + "\t".join(region))

    return output_table


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~CALCULATING CHANGES IN OUT DEGREE BY TF EDGES~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


def tf_edge_delta_out(
    crc_folder,
    bam_list,
    analysis_name,
    edge_table_path_1,
    edge_table_path_2,
    group1_list,
    group2_list,
    output="",
):
    """Calculates changes in group out degree at each predicted motif occurrence (by subpeaks)."""
    crc_folder = utils.format_folder(crc_folder, True)
    edge_path = merge_edge_tables(
        edge_table_path_1,
        edge_table_path_2,
        os.path.join(crc_folder, "{}_EDGE_TABLE.txt".format(analysis_name)),
    )

    # make a gff of the edge table
    edge_table = utils.parse_table(edge_path, "\t")
    edge_gff = []
    for line in edge_table[1:]:
        gff_line = [
            line[2],
            "{}_{}".format(line[0], line[1]),
            "",
            line[3],
            line[4],
            "",
            ".",
            "",
            "{}_{}".format(line[0], line[1]),
        ]
        edge_gff.append(gff_line)

    edge_gff_path = os.path.join(crc_folder, "{}_EDGE_TABLE.gff".format(analysis_name))
    utils.unparse_table(edge_gff, edge_gff_path, "\t")

    # direct the output to the crc folder
    signal_path = os.path.join(
        crc_folder, "{}_EDGE_TABLE_signal.txt".format(analysis_name)
    )

    all_group_list = group1_list + group2_list
    if not utils.check_output(signal_path, 0, 0):
        signal_table_list = pipeline_utils.map_regions(
            bam_list,
            [edge_gff_path],
            crc_folder,
            crc_folder,
            all_group_list,
            True,
            signal_path,
            extend_reads_to=100,
        )
        print(signal_table_list)
    else:
        print("Found previous signal table at {}".format(signal_path))

    # now bring in the signal table as a dictionary using the locus line as the id
    print("making log2 group1 vs group2 signal table at edges")
    signal_table = utils.parse_table(signal_path, "\t")

    # figure out columns for group1 and group2
    group1_columns = [signal_table[0].index(name) for name in group1_list]
    group2_columns = [signal_table[0].index(name) for name in group2_list]
    group1_signal_vector = []
    group2_signal_vector = []
    for line in signal_table[1:]:
        group1_signal = numpy.mean([float(line[col]) for col in group1_columns])
        group2_signal = numpy.mean([float(line[col]) for col in group2_columns])

        group1_signal_vector.append(group1_signal)
        group2_signal_vector.append(group2_signal)

    group1_median = numpy.median(group1_signal_vector)
    group2_median = numpy.median(group2_signal_vector)

    print("group1 median signal")
    print(group1_median)
    print("group2 median signal")
    print(group2_median)

    # now that we have the median, we can take edges where at least 1 edge is above the median
    # and both are above zero and generate a new table w/ the fold change
    signal_filtered_path = signal_path.replace(".txt", "_filtered.txt")
    if utils.check_output(signal_filtered_path, 0, 0):
        print(
            "Found filtered signal table for edges at {}".format(signal_filtered_path)
        )
        signal_table_filtered = utils.parse_table(signal_filtered_path, "\t")
    else:
        signal_table_filtered = [
            signal_table[0] + ["GROUP1_MEAN", "GROUP2_MEAN", "GROUP1_vs_GROUP2_LOG2"]
        ]
        for line in signal_table[1:]:
            group1_signal = numpy.mean([float(line[col]) for col in group1_columns])
            group2_signal = numpy.mean([float(line[col]) for col in group2_columns])

            if (group1_signal > group1_median or group2_signal > group2_median) and min(
                group1_signal, group2_signal
            ) > 0:
                delta = numpy.log2(group1_signal / group2_signal)
                new_line = line + [group1_signal, group2_signal, delta]
                signal_table_filtered.append(new_line)

        utils.unparse_table(signal_table_filtered, signal_filtered_path, "\t")

    # now get a list of all TFs in the system
    tf_list = utils.uniquify(
        [line[0].split("_")[0] for line in signal_table_filtered[1:]]
    )
    tf_list.sort()
    print(tf_list)

    out_degree_table = [
        [
            "TF_NAME",
            "EDGE_COUNT",
            "DELTA_MEAN",
            "DELTA_MEDIAN",
            "DELTA_STD",
            "DELTA_SEM",
        ]
    ]

    for tf_name in tf_list:
        print(tf_name)
        edge_vector = [
            float(line[-1])
            for line in signal_table_filtered[1:]
            if line[0].split("_")[0] == tf_name
        ]

        edge_count = len(edge_vector)
        delta_mean = round(numpy.mean(edge_vector), 4)
        delta_median = round(numpy.median(edge_vector), 4)
        delta_std = round(numpy.std(edge_vector), 4)
        delta_sem = round(stats.sem(edge_vector), 4)
        tf_out_line = [
            tf_name,
            edge_count,
            delta_mean,
            delta_median,
            delta_std,
            delta_sem,
        ]
        out_degree_table.append(tf_out_line)

    # set final output
    if not output:
        output_path = os.path.join(
            crc_folder, "{}_EDGE_DELTA_OUT.txt".format(analysis_name)
        )
    else:
        output_path = output

    utils.unparse_table(out_degree_table, output_path, "\t")
    print(output_path)
    return output_path


# ==========================================================================
# ===========================MAIN METHOD====================================
# ==========================================================================


def main():

    print("main analysis")
    args = parse_args()

    print("changing directory to project folder")
    os.chdir(args.output)

    analysis_name = args.name
    bams_1 = args.bams_1
    bams_2 = args.bams_2
    bams_list = bams_1 + bams_2
    edge_table_1 = args.edge_table_1
    edge_table_2 = args.edge_table_2
    group_1_list = args.names_1 or [
        os.path.basename(bam).split(".")[0] for bam in bams_1
    ]
    group_2_list = args.names_2 or [
        os.path.basename(bam).split(".")[0] for bam in bams_2
    ]

    print("\n\n")
    print("#======================================================================")
    print("#========================I. CHECKING SEQ DATA==========================")
    print("#======================================================================")
    print("\n\n")

    # Check that each edge table and all bam and .bai files are accessible
    for bam in bams_list:
        if not os.path.isfile(bam):
            raise FileNotFoundError("BAM file {} does not exist.".format(bam))
        if not os.path.isfile("{}.bai".format(bam)):
            raise FileNotFoundError(
                "BAM index (.bai) file {} does not exist.".format(bam)
            )
        print("{} and its index exist".format(bam))
    for table in [edge_table_1, edge_table_2]:
        if not os.path.isfile(table):
            raise FileNotFoundError("Edge table {} does not exist.".format(table))
        print("{} exists".format(table))
    print("\nEverything looks OK")

    print("\n\n")
    print("#======================================================================")
    print("#========================II. DELTA OUT BY EDGE=========================")
    print("#======================================================================")
    print("\n\n")

    tf_edge_delta_out(
        args.output,
        bams_list,
        analysis_name,
        edge_table_1,
        edge_table_2,
        group_1_list,
        group_2_list,
    )


if __name__ == "__main__":
    main()
