#!/usr/bin/python

"""
The MIT License (MIT)

Copyright (c) 2020 Charles Lin

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

# program to perform 2D clustering by enhancer signal
# can perform initial enhancer mapping or draw from a set of finished rose outputs

import argparse
import os
import sys
from collections import defaultdict

import numpy

from pipeline_tools.definitions import ROOT_DIR
from pipeline_tools.utils import pipeline_utils, utils

# ==========================================================================
# ==============================FUNCTIONS===================================
# ==========================================================================


def make_name_dict(data_file, rose_folder, names_list=[], enhancer_type="super"):
    """For each name, check for the presence of an enriched file or allEnhancer table.

    These are the files required for enhancer clustering.

    """
    data_dict = pipeline_utils.load_data_table(data_file)

    # draw the parent folder from the data_file
    parent_folder = utils.format_folder(
        os.path.dirname(os.path.abspath(data_file)), False
    )
    if parent_folder.count("data_tables") == 1:
        parent_folder = parent_folder.replace("data_tables/", "")
    print("Using {} as the parent folder".format(parent_folder))

    # check to see if a rose folder exists already
    if utils.format_folder(rose_folder, False):
        rose_exists = True
        rose_folder = utils.format_folder(rose_folder, False)
    else:
        rose_exists = False
        rose_folder = utils.format_folder(rose_folder, True)

    # check names_list to see if datasets exist
    if len(names_list) == 0:
        names_list = [
            name
            for name in data_dict
            if name.upper().count("WCE") == 0 and name.upper().count("INPUT") == 0
        ]
        # if no names_list is given, this filters out WCE

    # now check that all of the datasets at a minimum have a rose output OR enriched region file

    name_dict = defaultdict(dict)
    for name in names_list:
        # check if each dataset has a background
        background_name = data_dict[name]["background"]
        if background_name in data_dict:
            name_dict[name]["background"] = True
        else:
            name_dict[name]["background"] = False

        # assumes standard folder structure for enriched file
        enriched_file = os.path.join(
            parent_folder, "macsEnriched", data_dict[name]["enrichedMacs"]
        )

        print("Looking for macs output at {}".format(enriched_file))

        try:
            open(enriched_file, "r").close()
            name_dict[name]["enriched_file"] = enriched_file
        except (IOError, FileNotFoundError):
            name_dict[name]["enriched_file"] = ""

        # roseOutput looks for standard format rose output
        # need an allEnhancers table and a region table to proceed
        # if the rose folder doesn't exist, don't bother
        if rose_exists:
            try:
                rose_output_files = os.listdir(
                    os.path.join(rose_folder, "{}_ROSE".format(name))
                )
                if enhancer_type == "super":
                    enhancer_string = "AllEnhancers.table.txt"
                if enhancer_type == "stretch":
                    enhancer_string = "AllEnhancers_Length.table.txt"
                if enhancer_type == "superstretch":
                    enhancer_string = "AllEnhancers_SuperStretch.table.txt"

                all_enhancer_file_list = [
                    x
                    for x in rose_output_files
                    if x.count(enhancer_string) == 1 and x[0] != "."
                ]  # no weird hidden or temp files
                if all_enhancer_file_list:
                    name_dict[name]["enhancer_file"] = os.path.join(
                        rose_folder, "{}_ROSE".format(name), all_enhancer_file_list[0]
                    )
                else:
                    name_dict[name]["enhancer_file"] = ""
            except (OSError, FileNotFoundError):
                name_dict[name]["enhancer_file"] = ""
        else:
            name_dict[name]["enhancer_file"] = ""

        if (
            name_dict[name]["enhancer_file"] == ""
            and name_dict[name]["enriched_file"] == ""
        ):
            print(
                "INSUFFICIENT DATA TO RUN ENAHNCER ANALYSIS ON {}. PLEASE MAKE SURE ROSE OUTPUT "
                "OR MACS ENRICHED REGION PEAKS FILE EXISTS".format(name)
            )
            print(name_dict[name])
            sys.exit()

    return name_dict


def launch_enhancer_mapping(
    data_file,
    name_dict,
    output_folder,
    rose_folder,
    stitch,
    tss_distance,
    enhancer_type,
    mask_file="",
):
    """Launches enhancer mapping if needed from enriched region files."""
    names_list = list(name_dict.keys())

    # check to see if everything is good, if so return True and call it a day
    if len([x for x in names_list if len(name_dict[x]["enhancer_file"]) > 0]) == len(
        names_list
    ):
        print("ENHANCER FILE OUTPUT FOUND FOR ALL DATASETS")
        return name_dict

    # if not, have to call rose
    rose_output_folder = utils.format_folder(rose_folder, True)

    queue_list = []
    for name in names_list:
        # check to see if we need to call rose
        if name_dict[name]["enhancer_file"] == "":
            # get the enriched file
            enriched_file = name_dict[name]["enriched_file"]
            # call rose
            print("CALLING ROSE FOR {}".format(name))
            bash_file_name = pipeline_utils.call_rose2(
                data_file,
                "",
                rose_output_folder,
                [name],
                [],
                enriched_file,
                tss_distance,
                stitch,
                mask=mask_file,
            )
            print(bash_file_name)
            os.system("bash {}".format(bash_file_name))
            # add name to queue list
            queue_list.append(name)

    # define the enhancer type
    if enhancer_type == "super":
        enhancer_string = "AllEnhancers.table.txt"
    if enhancer_type == "stretch":
        enhancer_string = "AllEnhancers_Length.table.txt"
    if enhancer_type == "superstretch":
        enhancer_string = "AllEnhancers_SuperStretch.table.txt"

    # now check for completion of datasets
    for name in queue_list:
        # check for the AllEnhancers table
        enhancer_file = os.path.join(
            rose_output_folder,
            "{}_ROSE".format(name),
            "{}_peaks_{}".format(name, enhancer_string),
        )

        print("CHECKING FOR {} ROSE OUTPUT IN {}".format(name, enhancer_file))
        if utils.check_output(enhancer_file, 1, 10):

            print("FOUND ENHANCER OUTPUT FOR {}".format(name))
            name_dict[name]["enhancer_file"] = enhancer_file
        else:
            # try finding it w/ a different name
            # this will bug out if nothing is there
            rose_folder = os.path.join(rose_output_folder, "{}_ROSE".format(name))
            rose_file_list = [
                x for x in os.listdir(rose_folder) if x[0] != "."
            ]  # no hidden files
            if not rose_file_list:
                print("No files found in {}".format(rose_folder))
                sys.exit()
            enhancer_file = pipeline_utils.get_file(enhancer_string, rose_file_list, rose_folder)
            name_dict[name]["enhancer_file"] = enhancer_file

    return name_dict


def make_median_dict(name_dict):
    """For each dataset returns the median background subtracted enhancer signal."""
    median_dict = {}
    for name in name_dict:
        # open up the allenhancer_table
        enhancer_table = utils.parse_table(name_dict[name]["enhancer_file"], "\t")

        if name_dict[name]["background"] is True:
            # assume header ends after line 5
            enhancer_vector = [
                float(line[6]) - float(line[7]) for line in enhancer_table[6:]
            ]
        else:
            enhancer_vector = [float(line[6]) for line in enhancer_table[6:]]

        median_dict[name] = numpy.median(enhancer_vector)

    return median_dict


def make_se_collection(enhancer_file, name, super_only=True):
    """Return a locus collection from a super table."""
    enhancer_table = utils.parse_table(enhancer_file, "\t")
    enhancer_loci = []
    for line in enhancer_table:
        if line[0][0] == "#" or line[0][0] == "R":
            continue
        else:
            if super_only and int(line[-1]) == 0:
                break
            enhancer_loci.append(
                utils.Locus(
                    line[1], line[2], line[3], ".", "{}_{}".format(name, line[0])
                )
            )

    return utils.LocusCollection(enhancer_loci, 50)


def merge_collections(name_dict, analysis_name, output="", super_only=True):
    """Merge them collections."""
    all_loci = []
    names_list = list(name_dict.keys())
    for name in names_list:
        se_collection = make_se_collection(
            name_dict[name]["enhancer_file"], name, super_only
        )
        if super_only:
            print(
                "DATASET: {} HAS {} SUPERENHANCERS".format(
                    name, str(len(se_collection))
                )
            )
        else:
            print("DATASET: {} HAS {} ENHANCERS".format(name, str(len(se_collection))))
        all_loci += se_collection.get_loci()

    print(str(len(all_loci)))

    merged_collection = utils.LocusCollection(all_loci, 50)

    # stitch the collection together
    stitched_collection = merged_collection.stitch_collection()
    stitched_loci = list(stitched_collection.get_loci())
    print("IDENTIFIED {} CONSENSUS ENHANCER REGIONS".format(str(len(stitched_loci))))

    # sort by size and provide a unique ID
    size_list = [locus.len() for locus in stitched_loci]
    size_order = utils.order(size_list, decreasing=True)
    ordered_loci = [stitched_loci[i] for i in size_order]
    for i in range(len(ordered_loci)):
        ordered_loci[i].id = "merged_{}_{}".format(analysis_name, str(i + 1))

    merged_gff = []
    for locus in ordered_loci:
        new_line = [
            locus.chr,
            locus.id,
            "",
            locus.start,
            locus.end,
            "",
            locus.sense,
            "",
            locus.id,
        ]
        merged_gff.append(new_line)

    if len(output) == 0:
        return merged_gff
    else:
        print("writing merged gff to {}".format(output))
        utils.unparse_table(merged_gff, output, "\t")
        return output


def map_merged_gff(
    data_file, name_dict, merged_gff_file, analysis_name, output_folder, mask_file
):
    """Call rose on the merged_gff_file for all datasets."""
    data_dict = pipeline_utils.load_data_table(data_file)
    rose_parent_folder = os.path.join(output_folder, "rose")
    utils.format_folder(rose_parent_folder, True)
    gff_name = os.path.basename(merged_gff_file).split(".")[0]
    bash_file_name = os.path.join(
        output_folder, "rose", "{}_roseCall.sh".format(analysis_name)
    )
    # names_list is just the first dataset
    # extrmap will have to have all other datasets + their backgrounds

    names_list = list(name_dict.keys())
    names_list.sort()
    extra_map = []
    for name in names_list[1:]:
        if name_dict[name]["background"]:
            background_name = data_dict[name]["background"]
            if background_name in data_dict:
                extra_map += [name, background_name]
            else:
                print(
                    "ERROR: UNABLE TO FIND LISTED BACKGROUND DATASET {} FOR {}".format(
                        background_name, name
                    )
                )
                sys.exit()
        else:
            extra_map += [name]

    print(extra_map)

    # first check to see if this has already been done
    merged_region_map = os.path.join(
        output_folder,
        "rose",
        "{}_ROSE".format(names_list[0]),
        "{}_0KB_STITCHED_ENHANCER_REGION_MAP.txt".format(gff_name),
    )
    print("LOOKING FOR REGION MAP AT {}".format(merged_region_map))

    if utils.check_output(merged_region_map, 1, 1):
        print("FOUND PREVIOUS REGION MAP")

        return merged_region_map

    bash_file_name = pipeline_utils.call_rose2(
        data_file,
        "",
        rose_parent_folder,
        [names_list[0]],
        extra_map,
        merged_gff_file,
        0,
        0,
        bash_file_name,
        mask=mask_file,
    )

    bash_command = "bash {}".format(bash_file_name)
    os.system(bash_command)
    print("Running enhancer mapping command:\n{}".format(bash_command))

    if utils.check_output(merged_region_map, 5, 60):
        return merged_region_map
    else:
        print(
            "UNABLE TO CALL ROSE ENHANCER MAPPING ON CONSENSUS ENHANCER FILE {}.\nEXITING NOW"
            "".format(merged_gff_file)
        )
        sys.exit()


def make_enhancer_signal_table(
    name_dict, merged_region_map, median_dict, analysis_name, genome, output_folder
):
    """Makes a signal table.

    Each row is an enhancer and each column is the log2 background corrected signal vs. median.

    """
    # load in the region map
    region_map = utils.parse_table(merged_region_map, "\t")
    names_list = list(name_dict.keys())
    names_list.sort()
    signal_table = [
        ["REGION_ID", "CHROM", "START", "STOP", "NUM_LOCI", "CONSTITUENT_SIZE"]
        + names_list
    ]

    print("len of {} for names_list".format(len(names_list)))
    print(names_list)

    for line in region_map[1:]:
        new_line = line[0:6]
        # a little tricky here to add datasets sequentially
        i = 6  # start w/ the first column w/ data
        for name in names_list:
            if name_dict[name]["background"] is True:
                enhancer_index = int(i)
                i += 1
                control_index = int(i)
                i += 1
                try:
                    enhancer_signal = float(line[enhancer_index]) - float(
                        line[control_index]
                    )
                except IndexError:
                    print(line)
                    print(len(line))
                    print(enhancer_index)
                    print(control_index)
                    sys.exit()

            else:
                enhancer_index = int(i)
                i += 1
                enhancer_signal = float(line[enhancer_index])

            if enhancer_signal < 0:
                enhancer_signal = 0
            enhancer_signal = enhancer_signal / median_dict[name]
            new_line.append(enhancer_signal)

        signal_table.append(new_line)

    output_file = os.path.join(
        output_folder, "{}_{}_signal_table.txt".format(genome, analysis_name)
    )
    print("WRITING MEDIAN NORMALIZED SIGNAL TABLE TO {}".format(output_file))
    utils.unparse_table(signal_table, output_file, "\t")

    return output_file


def call_r_script(genome, output_folder, analysis_name, signal_table_file):
    """Call the R script to do clustering and heatmap."""
    cluster_table = os.path.join(
        output_folder, "{}_{}_clusterTable.txt".format(genome, analysis_name)
    )

    r_cmd = "Rscript {} {} {} {} {}".format(
        os.path.join(ROOT_DIR, "scripts", "clusterEnhancer.R"),
        genome,
        output_folder + '/',  # TODO: fix R script so it does not require '/'
        analysis_name,
        signal_table_file,
    )
    print("Calling command {}".format(r_cmd))

    os.system(r_cmd)

    print("Checking for cluster table output at {}".format(cluster_table))
    if utils.check_output(cluster_table, 1, 30):

        return cluster_table

    else:
        print("ERROR: CLUSTERING TABLE FAILED TO GENERATE")
        sys.exit()


# ==========================================================================
# =============================MAIN METHOD==================================
# ==========================================================================


def main():
    """Main function call."""
    parser = argparse.ArgumentParser()
    # required flags
    parser.add_argument(
        "-d",
        "--data",
        dest="data",
        default=None,
        required=True,
        help="Enter a data file for datasets to be processed",
    )
    parser.add_argument(
        "-o",
        "--output",
        dest="output",
        default=None,
        required=True,
        help="specify an output folder to write results to",
    )

    # additional args
    parser.add_argument(
        "-i",
        "--input",
        dest="input",
        required=False,
        help="Enter a comma separated list of names to analyze. Default will be all datasets",
    )

    parser.add_argument(
        "-n",
        "--name",
        dest="name",
        required=False,
        help="Enter a name for the analysis",
    )

    parser.add_argument(
        "-r",
        "--rose",
        dest="rose",
        required=False,
        help="Enter a folder to detect or write rose output",
    )

    parser.add_argument(
        "-a",
        "--all",
        dest="all",
        action="store_true",
        default=False,
        help="flag to run analysis on ALL enhancers (this is much slower)",
    )
    parser.add_argument(
        "-s",
        "--stitch",
        dest="stitch",
        default="",
        help=(
            "specify a fixed stitch distance for all datasets, otherwise will compute stitching "
            "automatically on each dataset"
        ),
    )
    parser.add_argument(
        "-e",
        "--enhancer-type",
        dest="enhancer_type",
        default="super",
        help="specify type of enhancer to analyze: super, stretch, superStretch",
    )

    parser.add_argument(
        "-t",
        "--tss",
        dest="tss",
        default=2500,
        help="specify a tss exclusion window. default is 2500bp",
    )

    parser.add_argument(
        "--mask",
        dest="mask",
        required=False,
        help="Create a mask set of regions to filter out of analysis. must be .bed or .gff format",
    )

    args = parser.parse_args()
    print(args)

    # pull in the data_file and create a data_dict
    data_file = args.data

    # now the output folder
    output_folder = utils.format_folder(
        args.output, True
    )  # check and create the output folder
    # now the rose folder
    if args.rose:
        rose_folder = args.rose
    else:
        rose_folder = os.path.join(output_folder, "rose")

    if args.input:
        names_list = args.input.split(",")
    else:
        names_list = []

    # get the genome
    data_dict = pipeline_utils.load_data_table(data_file)
    genome = data_dict[list(data_dict.keys())[0]]["genome"]

    # check if using only supers
    if args.all:
        super_only = False
    else:
        super_only = True

    # get the anlysis name
    if args.name:
        analysis_name = args.name
    else:
        analysis_name = "enhancers"

    # check for a stitching parameter
    if len(str(args.stitch)) > 0:
        stitch = str(args.stitch)
    else:
        stitch = ""

    # check for the tss parameter
    tss_distance = int(args.tss)

    # check enhancer type
    enhancer_type = args.enhancer_type.lower()
    if ["super", "superstretch", "stretch"].count(enhancer_type) == 0:
        print("ERROR: unsupported enhancer type {}".format(enhancer_type))
        sys.exit()

    # see if there's a mask
    if args.mask:
        mask_file = args.mask
    else:
        mask_file = ""

    # =====================================================
    # =================SUMMARIZE INPUTS====================
    # =====================================================

    print("WORKING IN GENOME {}".format(genome))
    print("DRAWING DATA FROM {} AND ROSE FOLDER {}".format(data_file, rose_folder))
    print("USING {} AS THE OUTPUT FOLDER".format(output_folder))

    # =====================================================
    # ==============ESTABLISH ALL WORKING FILES============
    # =====================================================

    print("\n\n\nESTABLISHING WORKING FILES")
    name_dict = make_name_dict(data_file, rose_folder, names_list, enhancer_type)

    print(name_dict)

    print("STARTING ANALYSIS ON THE FOLLOWING DATASETS:")
    print(list(name_dict.keys()))

    for name in name_dict:
        if len(name_dict[name]["enhancer_file"]) == 0:
            print("NO ROSE OUTPUT FOR {}".format(name))

    # =====================================================
    # ==============LAUNCH ENHANCER MAPPING================
    # =====================================================

    print("\n\n\nLAUNCHING ENHANCER MAPPING (IF NECESSARY)")
    name_dict = launch_enhancer_mapping(
        data_file,
        name_dict,
        output_folder,
        rose_folder,
        stitch,
        tss_distance,
        enhancer_type,
        mask_file,
    )
    print(name_dict)

    # =====================================================
    # ====================GET MEDIAN SIGNAL================
    # =====================================================

    print("\n\n\nGETTING MEDIAN ENHANCER SIGNAL FROM EACH SAMPLE")
    median_dict = make_median_dict(name_dict)

    print(median_dict)

    # =====================================================
    # ====================MERGING ENHANCERS================
    # =====================================================

    print("\n\n\nIDENTIFYING CONSENSUS ENHANCER REGIONS")

    merged_gff_file = os.path.join(
        output_folder, "{}_{}_-0_+0.gff".format(genome, analysis_name)
    )
    merged_gff_file = merge_collections(
        name_dict, analysis_name, merged_gff_file, super_only
    )

    # =====================================================
    # ===============MAP TO MERGED REGIONS=================
    # =====================================================

    print("\n\n\nMAPPING DATA TO CONSENSUS ENHANCER REGIONS")
    merged_region_map = map_merged_gff(
        data_file, name_dict, merged_gff_file, analysis_name, output_folder, mask_file
    )

    # =====================================================
    # ==============CORRECT FOR MEDIAN SIGNAL==============
    # =====================================================

    print("\n\n\nCREATING ENHANCER SIGNAL TABLE")
    signal_table_file = make_enhancer_signal_table(
        name_dict, merged_region_map, median_dict, analysis_name, genome, output_folder
    )

    # =====================================================
    # ===============CALL CLUSTERING R SCRIPT==============
    # =====================================================

    print("\n\n\nGENERATING CLUSTERING OUTPUT")
    cluster_table_file = call_r_script(
        genome, output_folder, analysis_name, signal_table_file
    )
    # output should be
    # png of cluster gram with rows as genes
    # png of cluster gram of samples w/ tree
    # ordered table w/ cluster assignment
    # similarity matrix for samples

    # =====================================================
    # =============GENE MAPPING BY CLUSTER=================
    # =====================================================

    cmd = "ROSE2_geneMapper -g {} -i {}".format(genome, cluster_table_file)
    os.system(cmd)

    print("FINISHED")


if __name__ == "__main__":
    main()
