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

# pipeline to run dynamic enhancer analysis

import argparse
import glob
import os
import sys
import time
from shutil import copyfile

import numpy

from pipeline_tools.definitions import ROOT_DIR
from pipeline_tools.utils import pipeline_utils, utils

# ================================================================================
# =================================FUNCTIONS======================================
# ================================================================================


def make_rose_dict(rose_folder):
    """Analyze a rose folder to try to find all of the various necessary files.

    Create a dictionary with their full paths.

    """
    print(rose_folder)

    if not utils.format_folder(rose_folder, False):
        print("Folder {} does not exist".format(rose_folder))
        sys.exit()

    rose_folder = utils.format_folder(rose_folder, False)
    rose_file_list = [
        x for x in os.listdir(rose_folder) if x[0] != "."
    ]  # no hidden files
    if not rose_file_list:
        print("No files found in {}".format(rose_folder))
        sys.exit()

    # create a dictionary to store stuff
    rose_dict = {}
    # there are 5 files that we're interested in
    # REGION_MAP, AllEnhancers.table.txt, SuperEnhancers.table.txt, ENHANCER_TO_GENE,
    # Enhancers_withSuper.bed

    # sequentially find each one and add the full path to the rose_dict
    rose_dict["AllEnhancer"] = pipeline_utils.get_file(
        "AllEnhancers.table.txt", rose_file_list, rose_folder
    )
    rose_dict["super"] = pipeline_utils.get_file(
        "SuperEnhancers.table.txt", rose_file_list, rose_folder
    )
    rose_dict["stretch"] = pipeline_utils.get_file(
        "_StretchEnhancers.table.txt", rose_file_list, rose_folder
    )
    rose_dict["superstretch"] = pipeline_utils.get_file(
        "SuperStretchEnhancers.table.txt", rose_file_list, rose_folder
    )

    rose_dict["enhancer_to_gene"] = pipeline_utils.get_file(
        "_SuperEnhancers_ENHANCER_TO_GENE", rose_file_list, rose_folder
    )
    rose_dict["RegionMap"] = pipeline_utils.get_file(
        "REGION_MAP", rose_file_list, rose_folder
    )
    rose_dict["bed"] = pipeline_utils.get_file(
        "Enhancers_withSuper.bed", rose_file_list, rose_folder
    )

    return rose_dict


def get_median_signal(enhancer_file, name, data_file):
    """Return the median enhancer signal of a file."""
    data_dict = pipeline_utils.load_data_table(data_file)
    enhancer_table = utils.parse_table(enhancer_file, "\t")
    background_name = data_dict[name]["background"]
    if background_name in data_dict:
        enhancer_vector = [
            float(line[6]) - float(line[7]) for line in enhancer_table[6:]
        ]
    else:
        enhancer_vector = [float(line[6]) for line in enhancer_table[6:]]

    median = numpy.median(enhancer_vector)

    return median


def make_se_collection(enhancer_file, name, top=0):
    """Return a locus collection from a super table.

    Top gives the number of rows.

    """
    enhancer_table = utils.parse_table(enhancer_file, "\t")
    super_loci = []

    ticker = 0
    for line in enhancer_table:
        if line[0][0] == "#" or line[0][0] == "R":
            continue
        else:
            ticker += 1
            super_loci.append(
                utils.Locus(
                    line[1], line[2], line[3], ".", "{}_{}".format(name, line[0])
                )
            )

            if ticker == top:
                break

    return utils.LocusCollection(super_loci, 50)


def make_se_dict(enhancer_file, name, super_only=True):
    """Make an attribute dict for enhancers keyed by uniqueID."""
    se_dict = {}
    enhancer_table = utils.parse_table(enhancer_file, "\t")
    for line in enhancer_table:
        if line[0][0] == "#":
            continue
        if line[0][0] == "R":
            header = line
            sup_column = header.index("isSuper")
            continue
        if super_only:
            if int(line[sup_column]) == 1:
                rank = int(line[-2])
                enhancer_id = "{}_{}".format(name, line[0])
                se_dict[enhancer_id] = {"rank": rank}
        else:
            rank = int(line[-2])
            enhancer_id = "{}_{}".format(name, line[0])
            se_dict[enhancer_id] = {"rank": rank}

    return se_dict


def merge_collections(super_file1, super_file2, name1, name2, output=""):
    """Merge them collections."""
    con_super_collection = make_se_collection(super_file1, name1)
    tnf_super_collection = make_se_collection(super_file2, name2)

    # now merge them
    merged_loci = con_super_collection.get_loci() + tnf_super_collection.get_loci()
    merged_collection = utils.LocusCollection(merged_loci, 50)

    # stitch the collection together
    stitched_collection = merged_collection.stitch_collection()
    stitched_loci = stitched_collection.get_loci()

    # loci that are in both get renamed with a new unique identifier
    renamed_loci = []
    ticker = 1
    for locus in stitched_loci:
        if len(con_super_collection.get_overlap(locus)) > 0 and len(
            tnf_super_collection.get_overlap(locus)
        ):
            new_id = "CONSERVED_{}".format(str(ticker))
            ticker += 1
            locus.id = new_id
        else:
            locus.id = locus.id[2:]
        renamed_loci.append(locus)

    # now we turn this into a gff and write it out
    gff = utils.locus_collection_to_gff(utils.LocusCollection(renamed_loci, 50))

    if len(output) == 0:
        return gff
    else:
        print("writing merged gff to {}".format(output))
        utils.unparse_table(gff, output, "\t")
        return output


def call_rose_merged(data_file, merged_gff_file, name1, name2, parent_folder):
    """Make a rose call for the merged supers."""
    data_dict = pipeline_utils.load_data_table(data_file)
    background_name1 = data_dict[name1]["background"]
    background_name2 = data_dict[name2]["background"]
    if background_name1 in data_dict and background_name2 in data_dict:
        has_background = True
    elif background_name1 not in data_dict and background_name2 not in data_dict:
        has_background = False
    else:
        print(
            "ERROR: Only 1 dataset has a background file. This is a very very bad idea"
        )
        sys.exit()

    if has_background:
        names_list = [name1]
        extra_map = [name2, data_dict[name2]["background"]]
    else:
        names_list = [name1]
        extra_map = [name2]

    return pipeline_utils.call_rose2(
        data_file,
        "",
        parent_folder,
        names_list,
        extra_map,
        merged_gff_file,
        tss=0,
        stitch=0,
    )


def call_merge_supers(
    data_file, super_file1, super_file2, name1, name2, merge_name, genome, parent_folder
):
    """Call ROSE2 on merged super enhancers."""
    merged_gff_file = "%s%s_%s_MERGED_REGIONS_-0_+0.gff" % (
        parent_folder,
        genome.upper(),
        merge_name,
    )

    # check to make sure this hasn't been done yet
    rose_output = os.path.join(
        parent_folder,
        "{}_ROSE".format(name1),
        "{}_{}_MERGED_REGIONS_-0_+0_SuperEnhancers_ENHANCER_TO_GENE.txt".format(
            genome.upper(), merge_name
        ),
    )

    try:
        utils.parse_table(rose_output, "\t")
        print("ROSE OUTPUT ALREADY FOUND HERE {}".format(rose_output))
        return rose_output
    except (FileNotFoundError, IOError):
        print(
            "MERGING ENHANCER REGIONS FROM {} and {}".format(super_file1, super_file2)
        )
        merged_gff = merge_collections(
            super_file1, super_file2, name1, name2, merged_gff_file
        )

        # call rose on the merged collection
        rose_bash_file = call_rose_merged(
            data_file, merged_gff, name1, name2, parent_folder
        )
        print(rose_bash_file)

        # run the bash command
        os.system("bash {}".format(rose_bash_file))

        # check for and return output
        if utils.check_output(rose_output, 1, 10):
            return rose_output
        else:
            # try finding it w/ a different name
            # this will bug out if nothing is there
            rose_folder = os.path.join(parent_folder, "{}_ROSE".format(name1))
            rose_file_list = [
                x for x in os.listdir(rose_folder) if x[0] != "."
            ]  # no hidden files
            if not rose_file_list:
                print("No files found in {}".format(rose_folder))
                sys.exit()

            pipeline_utils.get_file(
                "_SuperEnhancers_ENHANCER_TO_GENE.txt", rose_file_list, rose_folder
            )


def call_delta_r_script(
    merged_gff_file,
    parent_folder,
    data_file,
    name1,
    name2,
    all_file1,
    all_file2,
    median_scale=False,
):
    """Create R script command."""
    if median_scale:
        median1 = get_median_signal(all_file1, name1, data_file)
        median2 = get_median_signal(all_file2, name2, data_file)
        print("normalizing signal for {} by median value of {}".format(name1, median1))
        print("normalizing signal for {} by median value of {}".format(name2, median2))

    else:
        median1 = 1
        median2 = 1

    gff_name = os.path.basename(merged_gff_file).split(".")[0]
    stitched_file = os.path.join(
        parent_folder,
        "{}_ROSE".format(name1),
        "{}_0KB_STITCHED_ENHANCER_REGION_MAP.txt".format(gff_name),
    )

    rcmd = "Rscript {} {} {} {} {} {}".format(
        os.path.join(ROOT_DIR, "scripts", "dynamicEnhancer_plot.R"),
        stitched_file,
        name1,
        name2,
        median1,
        median2,
    )

    return rcmd


def call_rank_r_script(enhancer_rank_file, name1, name2, super_file1, super_file2):
    """Run the R script."""
    enhancer_collection1 = make_se_collection(super_file1, name1, False)
    enhancer_collection2 = make_se_collection(super_file2, name2, False)

    n_super1 = len(enhancer_collection1)
    n_super2 = len(enhancer_collection2)

    rcmd = "Rscript {} {} {} {} {} {}".format(
        os.path.join(ROOT_DIR, "scripts", "dynamicEnhancer_rank.R"),
        enhancer_rank_file,
        name1,
        name2,
        n_super1,
        n_super2,
    )

    return rcmd


def call_rose_gene_mapper(merged_gff_file, genome, parent_folder, name1):
    """Call the rose gene mapper w/ 100kb window."""
    gff_name = os.path.basename(merged_gff_file).split(".")[0]
    stitched_file = os.path.join(
        parent_folder,
        "{}_ROSE".format(name1),
        "{}_0KB_STITCHED_ENHANCER_REGION_MAP.txt".format(gff_name),
    )

    delta_file = stitched_file.replace("REGION_MAP", "DELTA")

    cmd = "ROSE2_geneMapper -g {} -i {} -w 100000".format(genome, delta_file)
    os.system(cmd)
    print(cmd)


def assign_enhancer_rank(
    enhancer_to_gene_file, enhancer_file1, enhancer_file2, name1, name2, rank_output=""
):
    """Assign enhancer rank to genes.

    For all genes in the enhancer_to_gene table, assign the highest overlapping ranked enhancer
    in the other tables.

    """
    enhancer_to_gene = utils.parse_table(enhancer_to_gene_file, "\t")

    enhancer_collection1 = make_se_collection(enhancer_file1, name1, False)
    enhancer_collection2 = make_se_collection(enhancer_file2, name2, False)

    enhancer_dict1 = make_se_dict(enhancer_file1, name1, False)
    enhancer_dict2 = make_se_dict(enhancer_file2, name2, False)

    # we're going to update the enhancer_to_gene_table
    enhancer_to_gene[0] += ["{}_rank".format(name1), "{}_rank".format(name2)]
    for i in range(1, len(enhancer_to_gene)):
        line = enhancer_to_gene[i]
        locus_line = utils.Locus(line[1], line[2], line[3], ".", line[0])

        # if the enhancer doesn't exist, its ranking is dead last on the enhancer list
        enhancer1_overlap = enhancer_collection1.get_overlap(locus_line, "both")
        if len(enhancer1_overlap) == 0:
            enhancer1_rank = len(enhancer_collection1)
        else:
            rank_list1 = [enhancer_dict1[x.id]["rank"] for x in enhancer1_overlap]
            enhancer1_rank = min(rank_list1)

        enhancer2_overlap = enhancer_collection2.get_overlap(locus_line, "both")
        if len(enhancer2_overlap) == 0:
            enhancer2_rank = len(enhancer_collection2)
        else:
            rank_list2 = [enhancer_dict2[x.id]["rank"] for x in enhancer2_overlap]
            enhancer2_rank = min(rank_list2)
        enhancer_to_gene[i] += [enhancer1_rank, enhancer2_rank]

    if len(rank_output) == 0:
        return enhancer_to_gene
    else:
        utils.unparse_table(enhancer_to_gene, rank_output, "\t")


def finish_rank_output(
    data_file,
    rank_output,
    genome,
    merge_folder,
    merge_name,
    name1,
    name2,
    cut_off=1.5,
    window=100000,
    super_only=True,
    plot_bam=True,
):
    """Finish rank output.

    Clean up the rank output table. Make a gff of all of the gained/lost supers beyond a certain
    cut_off w/ a window. Make a list of gained genes and lost genes. Make a bed of gained loss.

    """
    data_dict = pipeline_utils.load_data_table(data_file)
    # making sure window and cut_off are int/float
    cut_off = float(cut_off)
    window = int(window)
    genome = genome.upper()

    # make the output folder
    output_folder = utils.format_folder(os.path.join(merge_folder, "output"), True)

    # bring in the old rank table
    rank_enhancer_table = utils.parse_table(rank_output, "\t")

    # make a new formatted table
    header = rank_enhancer_table[0]
    header[-4] = "DELTA RANK"
    header[-3] = "IS_SUPER"
    formatted_rank_table = [header]

    # the gffs
    gained_gff = []
    lost_gff = []

    gained_window_gff = []
    lost_window_gff = []

    if super_only:
        enhancer_type = "SUPERS"
    else:
        enhancer_type = "ENHANCERS"

    # the beds
    if super_only:
        gained_track_header = (
            'track name="{} {} only SEs" description="{} super enhancers that are found only in '
            '{} vs {}" itemRGB=On color=255,0,0'.format(
                genome, name2, genome, name2, name1
            )
        )
        gained_bed = [[gained_track_header]]
        conserved_track_header = (
            'track name="{} {} and {} SEs" description="{} super enhancers that are found in both'
            ' {} vs {}" itemRGB=On color=0,0,0'.format(
                genome, name1, name2, genome, name1, name2
            )
        )
        conserved_bed = [[conserved_track_header]]

        lost_track_header = (
            'track name="{} {} only SEs" description="{} super enhancers that are found only in '
            '{} vs {}" itemRGB=On color=0,255,0'.format(
                genome, name1, genome, name1, name2
            )
        )
        lost_bed = [[lost_track_header]]
    else:
        gained_track_header = (
            'track name="{} {} only enhancers" description="{} enhancers that are found only in '
            '{} vs {}" itemRGB=On color=255,0,0'.format(
                genome, name2, genome, name2, name1
            )
        )
        gained_bed = [[gained_track_header]]
        conserved_track_header = (
            'track name="{} {} and {} enhancers" description="{} enhancers that are found in both'
            ' {} vs {}" itemRGB=On color=0,0,0'.format(
                genome, name1, name2, genome, name1, name2
            )
        )
        conserved_bed = [[conserved_track_header]]

        lost_track_header = (
            'track name="{} {} only enhancers" description="{} enhancers that are found only in '
            '{} vs {}" itemRGB=On color=0,255,0'.format(
                genome, name1, genome, name1, name2
            )
        )
        lost_bed = [[lost_track_header]]

    # the genes
    gene_table = [
        [
            "GENE",
            "ENHANCER_ID",
            "ENHANCER_CHROM",
            "ENHANCER_START",
            "ENHANCER_STOP",
            header[6],
            header[7],
            header[8],
            "STATUS",
        ]
    ]

    for line in rank_enhancer_table[1:]:
        # fixing the enhancer ID
        line[0] = line[0].replace("_lociStitched", "")
        formatted_rank_table.append(line)

        # getting the genes
        gene_list = []
        gene_list += line[9].split(",")
        gene_list += line[10].split(",")
        gene_list += line[11].split(",")
        gene_list = [x for x in gene_list if len(x) > 0]
        gene_list = utils.uniquify(gene_list)
        gene_string = ",".join(gene_list)

        bed_line = [line[1], line[2], line[3], line[0], line[-4]]

        # for gained
        if float(line[6]) > cut_off:
            gff_line = [
                line[1],
                line[0],
                "",
                line[2],
                line[3],
                "",
                ".",
                "",
                gene_string,
            ]
            gff_window_line = [
                line[1],
                line[0],
                "",
                int(line[2]) - window,
                int(line[3]) + window,
                "",
                ".",
                "",
                gene_string,
            ]
            gained_gff.append(gff_line)
            gained_window_gff.append(gff_window_line)
            gene_status = name2
            gained_bed.append(bed_line)
        # for lost
        elif float(line[6]) < (-1 * cut_off):
            gff_line = [
                line[1],
                line[0],
                "",
                line[2],
                line[3],
                "",
                ".",
                "",
                gene_string,
            ]
            gff_window_line = [
                line[1],
                line[0],
                "",
                int(line[2]) - window,
                int(line[3]) + window,
                "",
                ".",
                "",
                gene_string,
            ]
            lost_gff.append(gff_line)
            lost_window_gff.append(gff_window_line)
            gene_status = name1
            lost_bed.append(bed_line)
        # for conserved
        else:
            gene_status = "CONSERVED"
            conserved_bed.append(bed_line)

        # now fill in the gene Table
        for gene in gene_list:
            gene_table_line = [
                gene,
                line[0],
                line[1],
                line[2],
                line[3],
                line[6],
                line[7],
                line[8],
                gene_status,
            ]
            gene_table.append(gene_table_line)

    # concat the bed
    full_bed = gained_bed + conserved_bed + lost_bed

    # start writing the output
    # there's the two gffs, the bed,the formatted table, the gene table

    # formatted table
    formatted_filename = os.path.join(
        output_folder,
        "{}_{}_MERGED_{}_RANK_TABLE.txt".format(genome, merge_name, enhancer_type),
    )
    utils.unparse_table(formatted_rank_table, formatted_filename, "\t")

    # gffs
    gff_folder = utils.format_folder(output_folder + "gff/", True)
    gff_filename_gained = os.path.join(
        gff_folder,
        "{}_{}_{}_ONLY_{}_-0_+0.gff".format(
            genome, merge_name, name2.upper(), enhancer_type
        ),
    )
    gff_filename_window_gained = os.path.join(
        gff_folder,
        "{}_{}_{}_ONLY_{}_-{}KB_+{}KB.gff".format(
            genome,
            merge_name,
            name2.upper(),
            enhancer_type,
            str(window // 1000),
            str(window // 1000),
        ),
    )

    gff_filename_lost = os.path.join(
        gff_folder,
        "{}_{}_{}_ONLY_{}_-0_+0.gff".format(
            genome, merge_name, name1.upper(), enhancer_type
        ),
    )
    gff_filename_window_lost = os.path.join(
        gff_folder,
        "{}_{}_{}_ONLY_{}_-{}KB_+{}KB.gff".format(
            genome,
            merge_name,
            name1.upper(),
            enhancer_type,
            str(window // 1000),
            str(window // 1000),
        ),
    )

    utils.unparse_table(gained_gff, gff_filename_gained, "\t")
    utils.unparse_table(gained_window_gff, gff_filename_window_gained, "\t")

    utils.unparse_table(lost_gff, gff_filename_lost, "\t")
    utils.unparse_table(lost_window_gff, gff_filename_window_lost, "\t")

    # bed
    bed_filename = os.path.join(
        output_folder, "{}_{}_MERGED_{}.bed".format(genome, merge_name, enhancer_type)
    )
    utils.unparse_table(full_bed, bed_filename, "\t")

    # gene_table
    gene_filename = os.path.join(
        output_folder,
        "{}_{}_MERGED_{}_GENE_TABLE.txt".format(genome, merge_name, enhancer_type),
    )
    utils.unparse_table(gene_table, gene_filename, "\t")

    # finally, move all of the plots to the output folder
    copyfile(
        glob.glob(os.path.join(merge_folder, "{}_ROSE".format(name1), "*.pdf"))[0],
        os.path.join(
            output_folder,
            "{}_{}_MERGED_{}_DELTA.pdf".format(genome, merge_name, enhancer_type),
        ),
    )

    copyfile(
        glob.glob(
            os.path.join(merge_folder, "{}_ROSE".format(name1), "*RANK_PLOT.png")
        )[0],
        os.path.join(
            output_folder,
            "{}_{}_MERGED_{}_RANK_PLOT.png".format(genome, merge_name, enhancer_type),
        ),
    )

    # now execute the bamPlot_turbo commands
    if plot_bam:
        bam1 = data_dict[name1]["bam"]
        bam2 = data_dict[name2]["bam"]
        bam_string = "{} {}".format(bam1, bam2)
        name_string = "{} {}".format(name1, name2)
        color_string = "0,0,0:100,100,100"

        if len(gained_gff) > 0:
            # gained command
            plot_title = "{}_ONLY_SE".format(name2)
            cmd = (
                "bamPlot_turbo -g {} -b {} -i {} -o {} -n {} -c {} -t {} -r -y UNIFORM -p "
                "MULTIPLE".format(
                    genome,
                    bam_string,
                    gff_filename_gained,
                    output_folder,
                    name_string,
                    color_string,
                    plot_title,
                )
            )
            os.system(cmd)

            # gained window command
            plot_title = "{}_ONLY_SE_{}KB_WINDOW".format(name2, str(window // 1000))
            cmd = (
                "bamPlot_turbo -g {} -b {} -i {} -o {} -n {} -c {} -t {} -r -y UNIFORM -p "
                "MULTIPLE".format(
                    genome,
                    bam_string,
                    gff_filename_window_gained,
                    output_folder,
                    name_string,
                    color_string,
                    plot_title,
                )
            )
            os.system(cmd)

        if len(lost_gff) > 0:
            # lost command
            plot_title = "{}_ONLY_SE".format(name1)
            cmd = (
                "bamPlot_turbo -g {} -b {} -i {} -o {} -n {} -c {} -t {} -r -y UNIFORM -p "
                "MULTIPLE".format(
                    genome,
                    bam_string,
                    gff_filename_lost,
                    output_folder,
                    name_string,
                    color_string,
                    plot_title,
                )
            )
            os.system(cmd)

            # lost command
            plot_title = "{}_ONLY_SE_{}KB_WINDOW".format(name1, str(window // 1000))
            cmd = (
                "bamPlot_turbo -g {} -b {} -i {} -o {} -n {} -c {} -t {} -r -y UNIFORM -p "
                "MULTIPLE".format(
                    genome,
                    bam_string,
                    gff_filename_window_lost,
                    output_folder,
                    name_string,
                    color_string,
                    plot_title,
                )
            )
            os.system(cmd)

    return


# ================================================================================
# ===============================MAIN RUN=========================================
# ================================================================================


def main():
    """Main run function."""
    parser = argparse.ArgumentParser()
    # required flags
    parser.add_argument(
        "-g",
        "--genome",
        dest="genome",
        required=True,
        help="Enter the genome build (HG18,HG19,MM9,RN4,RN6) for the project",
    )
    parser.add_argument(
        "-d",
        "--data",
        dest="data",
        required=True,
        help="Enter the data file for the project",
    )
    parser.add_argument(
        "-r",
        "--rose",
        dest="rose",
        required=True,
        help="Enter a comma separated list of rose folder",
    )
    parser.add_argument(
        "-o",
        "--output",
        dest="output",
        required=True,
        help="Enter the output folder for the project",
    )
    parser.add_argument(
        "-n",
        "--names",
        dest="names",
        required=True,
        help="Enter a comma separated list of names to go with the datasets",
    )

    # additional args
    parser.add_argument(
        "-p",
        "--plot",
        dest="plot",
        action="store_true",
        default=False,
        help="If flagged, will plot differential regions",
    )
    parser.add_argument(
        "-a",
        "--all",
        dest="all",
        action="store_true",
        default=False,
        help="If flagged, will run analysis for all enhancers and not just supers.",
    )
    parser.add_argument(
        "-m",
        "--median",
        dest="median",
        action="store_true",
        default=False,
        help="If flagged, will use median enhancer scaling",
    )
    parser.add_argument(
        "-e",
        "--enhancer-type",
        dest="enhancer_type",
        default="super",
        help="specify type of enhancer to analyze: super, stretch, superStretch",
    )

    args = parser.parse_args()

    print(args)

    genome = args.genome.upper()
    data_file = args.data

    rose_folder_string = args.rose
    rose_folder1, rose_folder2 = rose_folder_string.split(",")
    parent_folder = utils.format_folder(args.output, True)

    name_string = args.names
    name1, name2 = name_string.split(",")

    merge_name = "{}_{}_merged".format(name1, name2)

    # option for median scaling
    median_scale = args.median

    plot_bam = args.plot
    if args.all:
        super_only = False
    else:
        super_only = True

    if super_only and plot_bam:
        print(
            "Running dynamic enhancer analysis on all super enhancers in {} and {} and plotting "
            "output to {}".format(
                name1, name2, parent_folder
            )
        )
    if super_only and not plot_bam:
        print(
            "Running dynamic enhancer analysis on all super enhancers in {} and {} and writing "
            "output to {}".format(
                name1, name2, parent_folder
            )
        )
    if not super_only and plot_bam:
        print(
            "Running dynamic enhancer analysis on all enhancers in {} and {} and plotting output "
            "to {}. WARNING: Plotting all differential enhancers could take a while".format(
                name1, name2, parent_folder
            )
        )
    if not super_only and not plot_bam:
        print(
            "Running dynamic enhancer analysis on all enhancers in {} and {} and writing output "
            "to {}.".format(
                name1, name2, parent_folder
            )
        )

    # part 1
    print("PART1: analyzing ROSE output from {} and {}".format(name1, name2))
    # start with the all enhancer tables from the initial rose calls

    rose_folder1 = utils.format_folder(rose_folder1, False)
    rose_folder2 = utils.format_folder(rose_folder2, False)

    rose_dict1 = make_rose_dict(rose_folder1)
    rose_dict2 = make_rose_dict(rose_folder2)

    # choosing the type of enhancer to analyze
    enhancer_call_type = args.enhancer_type.lower()
    if super_only:
        print("ANALYZING ENHANCER TYPE: {}".format(enhancer_call_type.upper()))

    super_file1 = rose_dict1[enhancer_call_type]
    super_file2 = rose_dict2[enhancer_call_type]

    all_file1 = rose_dict1["AllEnhancer"]
    all_file2 = rose_dict2["AllEnhancer"]

    print("\tMERGING ENHANCERS AND CALLING ROSE")
    if super_only:
        if len(super_file1) == 0:
            print(
                "ERROR: UNABLE TO FIND {} FILES IN {}".format(
                    enhancer_call_type, rose_folder1
                )
            )
            sys.exit()
        if len(super_file2) == 0:
            print(
                "ERROR: UNABLE TO FIND {} FILES IN {}".format(
                    enhancer_call_type, rose_folder2
                )
            )
            sys.exit()
        rose_output = call_merge_supers(
            data_file,
            super_file1,
            super_file2,
            name1,
            name2,
            merge_name,
            genome,
            parent_folder,
        )

    else:
        rose_output = call_merge_supers(
            data_file,
            all_file1,
            all_file2,
            name1,
            name2,
            merge_name,
            genome,
            parent_folder,
        )

    print("\tCALCULATING ENHANCER DELTA AND MAKING PLOTS")

    # part2 is the R script
    merged_gff_file = os.path.join(
        parent_folder, "{}_{}_MERGED_REGIONS_-0_+0.gff".format(genome, merge_name)
    )
    rcmd = call_delta_r_script(
        merged_gff_file,
        parent_folder,
        data_file,
        name1,
        name2,
        all_file1,
        all_file2,
        median_scale,
    )
    print(rcmd)
    os.system(rcmd)

    time.sleep(30)
    call_rose_gene_mapper(merged_gff_file, genome, parent_folder, name1)

    # rank the genes

    # part 3
    # rank the delta
    print("PART 3: assinging ranks to differential enhancers")
    print("\tASSIGNING SUPER RANK TO MERGED ENHANCERS")

    gff_name = "{}_{}_MERGED_REGIONS_-0_+0".format(genome, merge_name)
    enhancer_to_gene_file = os.path.join(
        parent_folder,
        "{}_ROSE".format(name1),
        "{}_0KB_STITCHED_ENHANCER_DELTA_ENHANCER_TO_GENE_100KB.txt".format(gff_name),
    )
    if utils.check_output(enhancer_to_gene_file):
        rank_output = os.path.join(
            parent_folder,
            "{}_ROSE".format(name1),
            "{}_0KB_STITCHED_ENHANCER_DELTA_ENHANCER_TO_GENE_100KB_RANK.txt".format(
                gff_name
            ),
        )
        assign_enhancer_rank(
            enhancer_to_gene_file, all_file1, all_file2, name1, name2, rank_output
        )
    else:
        print("ERROR: DELTA SCRIPT OR ROSE GENE MAPPER FAILED TO RUN")
        sys.exit()

    # make the rank plot
    print("MAKING RANK PLOTS")
    if utils.check_output(rank_output):
        rcmd = call_rank_r_script(rank_output, name1, name2, super_file1, super_file2)
        print(rcmd)
        os.system(rcmd)
    else:
        print("ERROR: RANK PLOT SCRIPT FAILED TO RUN")
        sys.exit()

    time.sleep(30)

    print("FINISHING OUTPUT")
    finish_rank_output(
        data_file,
        rank_output,
        genome,
        parent_folder,
        merge_name,
        name1,
        name2,
        1,
        100000,
        super_only,
        plot_bam,
    )


if __name__ == "__main__":
    main()
