#!/usr/bin/env python

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

# Functions require samtools to be callable with the command 'samtools' and
# bamliquidator to be callable with the command 'bamliquidator'

# ==========================================================================
# =======================DEPENDENCIES=======================================
# ==========================================================================
import argparse
import os
import shutil
import sys
import subprocess

from pipeline_tools.definitions import ROOT_DIR
from pipeline_tools.utils import utils

# as of now the number of bins to sample the space is hard wired
n_bins = 200

# script that takes in a list of bams, makes an intermediate table file,
# and then calls R to make the plot
# updated for batching with bamliquidator magic

# ==========================================================================
# ======================HELPER FUNCTIONS====================================
# ==========================================================================


def load_annot_file(genome, gene_list=[]):
    """Load in the annotation.

    Create a gene_dict and transcription collection.

    """
    annotation_folder = os.path.join(ROOT_DIR, "annotation")
    genome_dict = {
        "HG18": os.path.join(annotation_folder, "hg18_refseq.ucsc"),
        "MM9": os.path.join(annotation_folder, "mm9_refseq.ucsc"),
        "MM10": os.path.join(annotation_folder, "mm10_refseq.ucsc"),
        "HG19": os.path.join(annotation_folder, "hg19_refseq.ucsc"),
        "HG19_RIBO": os.path.join(annotation_folder, "hg19_refseq.ucsc"),
        "RN4": os.path.join(annotation_folder, "rn4_refseq.ucsc"),
        "RN6": os.path.join(annotation_folder, "rn6_refseq.ucsc"),
        "HG38": os.path.join(annotation_folder, "hg38_refseq.ucsc"),
    }

    annot_file = genome_dict[genome.upper()]

    gene_dict = utils.make_genes(annot_file, gene_list, True)
    tx_collection = utils.make_transcript_collection(annot_file, 0, 0, 500, gene_list)

    return gene_dict, tx_collection


def taste_the_rainbow(n):
    """Samples rainbow color space."""
    from colorsys import hsv_to_rgb

    color_list = []
    n_range = [x // n for x in range(0, n)]
    for i in n_range:
        color = [int(255 * x) for x in list(hsv_to_rgb(i, 0.9, 0.9))]

        color_list.append(color)
    return color_list


def map_gff_line_to_annot(
    gff_line, out_folder, n_bins, gene_dict, tx_collection, sense="both", header=""
):
    """For every line produces a file with all of the rectangles to draw."""
    if not header:
        gff_string = "{}_{}_{}_{}".format(
            gff_line[0], gff_line[6], gff_line[3], gff_line[4]
        )
    else:
        gff_string = header
    diagram_table = [[0, 0, 0, 0]]
    name_table = [["", 0, 0]]
    gff_locus = utils.Locus(
        gff_line[0], int(gff_line[3]), int(gff_line[4]), gff_line[6], gff_line[1],
    )
    scale_factor = n_bins / gff_locus.len()
    # plotting buffer for diagrams
    plot_buffer = int(gff_locus.len() / n_bins * 20)

    overlap_loci = tx_collection.get_overlap(gff_locus, sense="both")
    gene_list = [locus.id for locus in overlap_loci]

    if gff_line[6] == "-":
        ref_point = int(gff_line[4])
    else:
        ref_point = int(gff_line[3])
    offset_collection = utils.LocusCollection([], 500)
    for gene_id in gene_list:

        gene = gene_dict[gene_id]

        print(gene.common_name())
        if len(gene.common_name()) > 1:
            name = gene.common_name()
        else:
            name = gene_id
        offset = 4 * len(offset_collection.get_overlap(gene.tx_locus()))
        offset_collection.append(
            utils.make_search_locus(gene.tx_locus(), plot_buffer, plot_buffer,)
        )
        # write the name of the gene down
        if gene.sense() == "+":
            gene_start = gene.tx_locus().start
        else:
            gene_start = gene.tx_locus().end
        gene_start = abs(gene_start - ref_point) * scale_factor
        name_table.append([name, gene_start, -2 - offset])

        # draw a line across the entire txLocus
        [start, stop] = [
            abs(x - ref_point) * scale_factor for x in gene.tx_locus().coords()
        ]
        diagram_table.append([start, -0.01 - offset, stop, 0.01 - offset])

        # now draw thin boxes for all tx_exons
        if gene.tx_exons():
            for tx_exon in gene.tx_exons():

                [start, stop] = [
                    abs(x - ref_point) * scale_factor for x in tx_exon.coords()
                ]

                diagram_table.append([start, -0.5 - offset, stop, 0.5 - offset])

        # now draw fatty boxes for the coding exons if any
        if gene.cd_exons():
            for cd_exon in gene.cd_exons():

                [start, stop] = [
                    abs(x - ref_point) * scale_factor for x in cd_exon.coords()
                ]

                diagram_table.append([start, -1 - offset, stop, 1 - offset])

    utils.unparse_table(
        diagram_table,
        os.path.join(out_folder, "{}_diagramTemp.txt".format(gff_string)),
        "\t",
    )
    utils.unparse_table(
        name_table,
        os.path.join(out_folder, "{}_nameTemp.txt".format(gff_string)),
        "\t",
    )


def make_bed_collection(bed_file_list):
    """Takes in a list of bed files and makes a single huge collection.

    Each locus has as its ID the name of the bed file.

    """
    bed_loci = []
    print("MAKING BED COLLECTION FOR:")
    for bed_file in bed_file_list:

        bed_name = os.path.basename(bed_file).split(".")[0]
        print(bed_name)
        bed = utils.parse_table(bed_file, "\t")
        for line in bed:
            if len(line) >= 3:
                # check that line[0]
                if line[0][0:3] == "chr":
                    try:
                        coords = [int(line[1]), int(line[2])]
                        bed_locus = utils.Locus(
                            line[0], min(coords), max(coords), ".", bed_name
                        )
                        bed_loci.append(bed_locus)
                    except ValueError:
                        pass

        print("IDENTIFIED {} BED REGIONS".format(str(len(bed_loci))))

    return utils.LocusCollection(bed_loci, 50)


def map_gff_line_to_bed(gff_line, out_folder, n_bins, bed_collection, header=""):
    """For every line produces a file with all of the rectangles to draw."""
    if not header:
        gff_string = "{}_{}_{}_{}".format(
            gff_line[0], gff_line[6], gff_line[3], gff_line[4]
        )
    else:
        gff_string = header
    diagram_table = [[0, 0, 0, 0]]
    name_table = [["", 0, 0]]
    gff_locus = utils.Locus(
        gff_line[0], int(gff_line[3]), int(gff_line[4]), gff_line[6], gff_line[1],
    )

    scale_factor = n_bins / gff_locus.len()

    overlap_loci = bed_collection.get_overlap(gff_locus, sense="both")
    print(
        "IDENTIFIED {} OVERLAPPING BED LOCI FOR REGION {}".format(
            str(len(overlap_loci)), gff_line,
        )
    )

    # since beds come from multiple sources, we want to figure out how to offset them
    offset_dict = {}  # this will store each ID name
    bed_names_list = utils.uniquify([locus.id for locus in overlap_loci])
    bed_names_list.sort()
    for i in range(len(bed_names_list)):
        offset_dict[bed_names_list[i]] = (
            2 * i
        )  # offsets different categories of bed regions

    if gff_line[6] == "-":
        ref_point = int(gff_line[4])
    else:
        ref_point = int(gff_line[3])

    # fill out the name table
    for name in bed_names_list:
        offset = offset_dict[name]
        name_table.append([name, 0, 0.0 - offset])

    for bed_locus in overlap_loci:

        offset = offset_dict[bed_locus.id]

        [start, stop] = [abs(x - ref_point) * scale_factor for x in bed_locus.coords()]

        diagram_table.append([start, -0.5 - offset, stop, 0.5 - offset])

    utils.unparse_table(
        diagram_table,
        os.path.join(out_folder, "{}_bedDiagramTemp.txt".format(gff_string)),
        "\t",
    )
    utils.unparse_table(
        name_table,
        os.path.join(out_folder, "{}_bedNameTemp.txt".format(gff_string)),
        "\t",
    )


def map_bam_to_gff_line(
    bam_file, mmr, name, gff_line, color, n_bins, sense="both", extension=200
):
    """Maps reads from a bam to a gff."""

    print("using a MMR/scaling denominator value of {}".format(mmr))

    line = gff_line[0:9]
    gff_locus = utils.Locus(line[0], int(line[3]), int(line[4]), line[6], line[1])

    # setting up the output clusterline
    color_line = color
    bam_name = os.path.basename(bam_file)
    cluster_line = [bam_name, gff_locus.id, name, gff_locus.__str__()] + color_line

    bin_size = gff_locus.len() // n_bins
    # some regions will be too short to get info on
    # we just kick these back and abandon them
    if not bin_size:
        cluster_line += ["NA"] * int(n_bins)
        return cluster_line

    # flippy flip if sense is negative
    sense_trans = str.maketrans("-+.", "+-+")
    if sense == "-":
        bam_sense = gff_locus.sense.translate(sense_trans)
    elif sense == "+":
        bam_sense = gff_locus.sense
    else:
        bam_sense = "."

    # using the bamliquidator to get the read_string
    bam_command = "bamliquidator {} {} {} {} {} {} {}".format(
        bam_file,
        gff_locus.chr,
        gff_locus.start,
        gff_locus.end,
        bam_sense,
        n_bins,
        extension,
    )
    get_reads = subprocess.Popen(
        bam_command,
        stdin=subprocess.PIPE,
        stderr=subprocess.PIPE,
        stdout=subprocess.PIPE,
        shell=True,
    )
    read_string = get_reads.communicate()
    den_list = read_string[0].decode("utf-8").split("\n")[:-1]

    # flip the denList if the actual gff region is -
    if gff_locus.sense == "-":
        den_list = den_list[::-1]

    # converting from units of total bp of read sequence per bin to rpm/bp
    den_list = [round(float(x) / bin_size / mmr, 4) for x in den_list]

    cluster_line += den_list

    return cluster_line


def call_r_plot(summary_file, out_file, y_scale, plot_style, multi):
    """Calls the R plotting thingy."""
    if multi is True:
        page_flag = "MULTIPLE_PAGE"
    else:
        page_flag = "SINGLE_PAGE"

    cmd = "Rscript {} {} {} {} {} {}".format(
        os.path.join(ROOT_DIR, "scripts", "bamPlot_turbo.R"),
        summary_file,
        out_file,
        y_scale,
        plot_style,
        page_flag,
    )
    print("calling command {}".format(cmd))
    return cmd


def make_bam_plot_tables(
    gff,
    genome,
    bam_file_list,
    color_list,
    n_bins,
    sense,
    extension,
    rpm,
    out_folder,
    names,
    title,
    bed_collection,
    scale=None,
):
    """Makes a plot table for each line of the gff mapped against all the bams in the bamList."""
    # load in the gff
    if isinstance(gff, str):
        gff = utils.parse_table(gff, "\t")

    # load in the annotation
    print("loading in annotation for {}".format(genome))
    gene_dict, tx_collection = load_annot_file(genome)

    # make an MMR dict so MMRs are only computed once
    print("Getting information about read depth in bams")
    mmr_dict = {}

    if scale:
        print("Applying scaling factors")
        scale_list = [float(x) for x in scale]
    else:
        scale_list = [1] * len(bam_file_list)

    # now iterate through the bam files
    for i, bam_file in enumerate(bam_file_list):
        # millionMappedReads
        idx_cmd = "samtools idxstats {}".format(bam_file)

        idx_pipe = subprocess.Popen(
            idx_cmd,
            stdin=subprocess.PIPE,
            stderr=subprocess.PIPE,
            stdout=subprocess.PIPE,
            shell=True,
        )  # TODO: this does not produce an error if samtools are not installed
        idx_stats = idx_pipe.communicate()
        idx_stats = idx_stats[0].decode("utf-8").split("\n")
        idx_stats = [line.split("\t") for line in idx_stats]
        raw_count = sum([int(line[2]) for line in idx_stats[:-1]])

        # implement scaling
        read_scale_factor = scale_list[i]

        if rpm:
            mmr = round(raw_count / 1000000 / read_scale_factor, 4)
        else:
            mmr = round(1 / read_scale_factor, 4)
        mmr_dict[bam_file] = mmr

    ticker = 1
    # go line by line in the gff
    summary_table = [
        [
            "DIAGRAM_TABLE",
            "NAME_TABLE",
            "BED_DIAGRAM_TABLE",
            "BED_NAME_TABLE",
            "PLOT_TABLE",
            "CHROM",
            "ID",
            "SENSE",
            "START",
            "END",
        ]
    ]
    for gff_line in gff:
        gff_string = "line_{}_{}_{}_{}_{}_{}".format(
            ticker, gff_line[0], gff_line[1], gff_line[6], gff_line[3], gff_line[4],
        )
        ticker += 1
        print("writing the gene diagram table for region {}".format(gff_line[1]))
        map_gff_line_to_annot(
            gff_line,
            out_folder,
            n_bins,
            gene_dict,
            tx_collection,
            sense="both",
            header=gff_string,
        )
        map_gff_line_to_bed(
            gff_line, out_folder, n_bins, bed_collection, header=gff_string,
        )
        out_table = []

        out_table.append(
            ["BAM", "GENE_ID", "NAME", "LOCUSLINE", "COLOR1", "COLOR2", "COLOR3"]
            + ["bin_" + str(n) for n in range(1, int(n_bins) + 1, 1)]
        )

        for i, bam_file in enumerate(bam_file_list):
            name = names[i]
            color = color_list[i]
            print(
                "getting data for location {} in dataset {}".format(
                    gff_line[1], bam_file
                )
            )
            mmr = mmr_dict[bam_file]
            new_line = map_bam_to_gff_line(
                bam_file, mmr, name, gff_line, color, n_bins, sense, extension,
            )
            out_table.append(new_line)

        # get the gene name
        if gff_line[1] in gene_dict:
            gene_name = gene_dict[gff_line[1]].common_name()
        else:
            gene_name = gff_line[1]
        utils.unparse_table(
            out_table,
            os.path.join(out_folder, "{}_plotTemp.txt".format(gff_string)),
            "\t",
        )
        diagram_table = os.path.join(
            out_folder, "{}_diagramTemp.txt".format(gff_string)
        )
        plot_table = os.path.join(out_folder, "{}_plotTemp.txt".format(gff_string))
        name_table = os.path.join(out_folder, "{}_nameTemp.txt".format(gff_string))
        bed_name_table = os.path.join(
            out_folder, "{}_bedNameTemp.txt".format(gff_string)
        )
        bed_diagram_table = os.path.join(
            out_folder, "{}_bedDiagramTemp.txt".format(gff_string)
        )
        summary_table.append(
            [
                diagram_table,
                name_table,
                bed_diagram_table,
                bed_name_table,
                plot_table,
                gff_line[0],
                gene_name,
                gff_line[6],
                gff_line[3],
                gff_line[4],
            ]
        )
    summary_table_file_name = os.path.join(out_folder, "{}_summary.txt".format(title))
    utils.unparse_table(summary_table, summary_table_file_name, "\t")
    return summary_table_file_name


def main():
    """Main run function."""
    parser = argparse.ArgumentParser()

    # required flags
    parser.add_argument(
        "-b",
        "--bam",
        dest="bam",
        nargs="*",
        help="Enter a comma/space separated list of .bam files to be processed.",
        required=True,
    )
    parser.add_argument(
        "-i",
        "--input",
        dest="input",
        type=str,
        help="Enter .gff or genomic region e.g. chr1:+:1-1000.",
        required=True,
    )
    parser.add_argument(
        "-g",
        "--genome",
        dest="genome",
        type=str,
        help="specify a genome, HG18,HG19,MM8,MM9,MM10 are currently supported",
        required=True,
    )

    # output flag
    parser.add_argument(
        "-o",
        "--output",
        dest="output",
        type=str,
        help="Enter the output folder.",
        required=True,
    )

    # additional options
    parser.add_argument(
        "--stretch-input",
        dest="stretch_input",
        default=None,
        type=int,
        help=(
            "Stretch the input regions to a minimum length in bp, e.g. 10000 (for"
            " 10kb)"
        ),
    )
    parser.add_argument(
        "-c",
        "--color",
        dest="color",
        default=None,
        nargs="*",
        help=(
            "Enter a colon or space separated list of colors e.g. "
            "255,0,0:255,125,0, default samples the rainbow"
        ),
    )
    parser.add_argument(
        "-s",
        "--sense",
        dest="sense",
        default="both",
        help="Map to '+','-' or 'both' strands. Default maps to both.",
    )
    parser.add_argument(
        "-e",
        "--extension",
        dest="extension",
        default=200,
        help="Extends reads by n bp. Default value is 200bp",
    )
    parser.add_argument(
        "-r",
        "--rpm",
        dest="rpm",
        action="store_true",
        default=False,
        help="Normalizes density to reads per million (rpm) Default is False",
    )
    parser.add_argument(
        "-y",
        "--yScale",
        dest="y_scale",
        default="relative",
        help=(
            "Choose either relative or uniform y axis scaling. options = "
            "'relative,uniform' Default is relative scaling"
        ),
    )
    parser.add_argument(
        "-n",
        "--names",
        dest="names",
        default=None,
        nargs="*",
        help="Enter a comma or space separated list of names for your bams",
    )
    parser.add_argument(
        "-p",
        "--plot",
        dest="plot",
        default="MULTIPLE",
        help=(
            "Choose either all lines on a single plot or multiple plots. options "
            "= 'SINGLE,MULTIPLE,MERGE'"
        ),
    )
    parser.add_argument(
        "-t",
        "--title",
        dest="title",
        default="",
        help=(
            "Specify a title for the output plot(s), default will be the "
            "coordinate region"
        ),
    )
    parser.add_argument(
        "-q",
        "--skip-cache",
        dest="skip_cache",
        action="store_true",
        default=False,
        help="Toggles option to skip loading annotation cache file",
    )

    parser.add_argument(
        "--scale",
        dest="scale",
        default=None,
        nargs="*",
        help=(
            "Enter a comma or space separated list of scaling factors for your "
            "bams. Default is none"
        ),
    )
    parser.add_argument(
        "--bed",
        dest="bed",
        nargs="*",
        help="Add a comma-delimited or space-delimited list of bed files to plot",
    )
    parser.add_argument(
        "--multi-page",
        dest="multi",
        action="store_true",
        default=False,
        help="If flagged will create a new pdf for each region",
    )

    # DEBUG OPTION TO SAVE TEMP FILES
    parser.add_argument(
        "--save-temp",
        dest="save",
        action="store_true",
        default=False,
        help="If flagged will save temporary files made by bamPlot",
    )

    args = parser.parse_args()

    print(args)

    if args.bam and args.input and args.genome and args.output:

        # Support a legacy mode where a ',' delimited multiple files
        bam_file_list = args.bam
        if len(bam_file_list) == 1:
            bam_file_list = bam_file_list[0].split(",")

        # Make sure these are actually files & readable (!)
        for filename in bam_file_list:
            assert os.access(filename, os.R_OK)

        # bringing in any beds
        if args.bed:
            bed_file_list = args.bed
            if len(bed_file_list) == 1:
                bed_file_list = bed_file_list[0].split(",")
            print(bed_file_list)
            bed_collection = make_bed_collection(bed_file_list)
        else:
            bed_collection = utils.LocusCollection([], 50)

        # Load the input for graphing. One of:
        # - A .gff
        # - A .bed
        # - a specific input region (e.g. chr10:.:93150000-93180000)

        valid_sense_options = {"+", "-", "."}
        if os.access(args.input, os.R_OK):
            if args.input.endswith(".bed"):
                # Uniquely graph every input of this bed
                parsed_input_bed = utils.parse_table(args.input, "\t")
                gff_name = os.path.basename(args.input)  # Graph title
                gff = None
                try:
                    if parsed_input_bed[0][5] in valid_sense_options:
                        # This .bed might have a sense parameter
                        gff = [
                            [e[0], "", args.input, e[1], e[2], "", e[5], "", ""]
                            for e in parsed_input_bed
                        ]
                except IndexError:
                    pass

                if gff is None:
                    print(
                        "Your bed doesn't have a valid sense parameter. Defaulting to both "
                        "strands, '.'"
                    )
                    # We only take chr/start/stop and ignore everything else.
                    gff = [
                        [e[0], "", args.input, e[1], e[2], "", ".", "", ""]
                        for e in parsed_input_bed
                    ]
            else:
                # Default to .gff, since that's the original behavior
                gff = utils.parse_table(args.input, "\t")
                gff_name = os.path.basename(args.input).split(".")[0]
        else:
            # means a coordinate line has been given e.g. chr1:+:1-100
            chrom_line = args.input.split(":")
            try:
                chrom = chrom_line[0]
                sense = chrom_line[1]
            except IndexError:
                print("Invalid input line or inaccessible file. Try: chr1:.:1-5000")
                exit()
            assert sense in valid_sense_options
            [start, end] = chrom_line[2].split("-")
            if chrom[0:3] != "chr":
                print("ERROR: UNRECOGNIZED GFF OR CHROMOSOME LINE INPUT")
                exit()
            gff_line = [chrom, "", args.input, start, end, "", sense, "", ""]
            gff_name = "{}_{}_{}_{}".format(chrom, sense, start, end)
            gff = [gff_line]

        # Consider stretching the regions to a fixed minimum size
        if args.stretch_input:
            print(
                "Stretching inputs to a minimum of: {} bp".format(
                    str(args.stretch_input)
                )
            )
            min_length = args.stretch_input
            stretch_gff = []
            for e in gff:
                difference = int(e[4]) - int(e[3])
                if difference < min_length:
                    pad = int((min_length - difference) / 2)
                    stretch_gff.append(
                        [
                            e[0],
                            e[1],
                            e[2],
                            int(e[3]) - pad,
                            int(e[4]) + pad,
                            e[5],
                            e[6],
                            e[7],
                            e[8],
                        ]
                    )
                else:
                    stretch_gff.append(e)

            gff = stretch_gff

        # Sanity test the gff object
        assert all([e[6] in valid_sense_options for e in gff])  # All strands are sane

        # bring in the genome
        genome = args.genome.upper()
        if not ["HG18", "HG19", "HG19_RIBO", "HG38", "MM9", "MM10", "RN4", "RN6"].count(
            genome
        ):
            print(
                "ERROR: UNSUPPORTED GENOME TYPE {}. USE HG19,HG18, RN4, MM9, or MM10".format(
                    genome,
                )
            )
            parser.print_help()
            exit()

        # bring in the rest of the options

        # output
        root_folder = args.output
        try:
            os.listdir(root_folder)
        except OSError:
            print("ERROR: UNABLE TO FIND OUTPUT DIRECTORY {}".format(root_folder))
            exit()

        # Get analysis title
        if not args.title:
            title = gff_name
        else:
            title = args.title

        # make a temp folder
        temp_folder = os.path.join(root_folder, title)
        print("CREATING TEMP FOLDER {}".format(temp_folder))
        utils.format_folder(temp_folder, create=True)

        # colors
        if args.color:
            color_list = args.color
            if len(color_list) == 1:
                color_list = color_list[0].split(":")
            color_list = [x.split(",") for x in color_list]
            if len(color_list) < len(bam_file_list):
                print(
                    "WARNING: FEWER COLORS THAN BAMS SPECIFIED. COLORS WILL BE RECYCLED"
                )
                # recycling the color list
                color_list += color_list * (len(bam_file_list) // len(color_list))
                color_list = color_list[: len(bam_file_list)]

        else:
            # cycles through the colors of the rainbow
            color_list = taste_the_rainbow(len(bam_file_list))

        # sense
        sense = args.sense

        extension = int(args.extension)

        rpm = args.rpm

        scale = args.scale
        if scale:
            if len(scale) == 1:
                scale = scale[0].split(",")

        y_scale = args.y_scale.upper()

        # names
        if args.names:
            names = args.names
            if len(names) == 1:
                names = names[0].split(",")

            if len(names) != len(bam_file_list):
                print("ERROR: NUMBER OF NAMES AND NUMBER OF BAMS DO NOT CORRESPOND")
                parser.print_help()
                exit()
        else:
            names = [os.path.basename(x) for x in bam_file_list]

        # plot style
        plot_style = args.plot.upper()
        if not ["SINGLE", "MULTIPLE", "MERGE"].count(plot_style):
            print("ERROR: PLOT STYLE {} NOT AN OPTION".format(plot_style))
            parser.print_help()
            exit()

        # now run!
        summary_table_file_name = make_bam_plot_tables(
            gff,
            genome,
            bam_file_list,
            color_list,
            n_bins,
            sense,
            extension,
            rpm,
            temp_folder,
            names,
            title,
            bed_collection,
            scale,
        )
        print("{} is the summary table".format(summary_table_file_name))

        # running the R command to plot
        multi = args.multi
        out_file = os.path.join(root_folder, "{}_plots.pdf".format(title))
        r_cmd = call_r_plot(
            summary_table_file_name, out_file, y_scale, plot_style, multi
        )

        # open a bash file
        bash_file_name = os.path.join(temp_folder, "{}_Rcmd.sh".format(title))
        with open(bash_file_name, "w") as bash_file:
            bash_file.write("#!/usr/bin/bash\n")
            bash_file.write(r_cmd)
        print("Wrote R command to {}".format(bash_file_name))
        os.system("bash {}".format(bash_file_name))

        # delete temp files
        if not args.save:
            if utils.check_output(out_file, 1, 10):
                # This is super dangerous (!). Add some sanity checks.
                assert " " not in temp_folder
                assert temp_folder != "/"
                shutil.rmtree(temp_folder)
                print("Removing temp folder: {}".format(temp_folder))
            else:
                print("ERROR: NO OUTPUT FILE {} DETECTED".format(out_file))

    else:
        parser.print_help()
        sys.exit()


if __name__ == "__main__":
    main()
