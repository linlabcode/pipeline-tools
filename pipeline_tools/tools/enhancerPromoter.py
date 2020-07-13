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

# Functions require bamliquidator_batch to be callable with a 'bamliquidator_batch' command

# ================================================================================
# =============================DEPENDENCIES=======================================
# ================================================================================

import argparse
import os
import re
import sys
from collections import defaultdict

import numpy
from pipeline_tools.definitions import ROOT_DIR
from pipeline_tools.utils import utils

# ================================================================================
# ============================GLOBAL PARAMETERS===================================
# ================================================================================

# using a paramater dictionary in liue of a yaml or json for now

param_dict = {
    "cpg_path": os.path.join(ROOT_DIR, "annotation", "hg19_cpg_islands.bed"),
}

# ================================================================================
# =================================FUNCTIONS======================================
# ================================================================================


def load_annot_file(genome, tss_window, gene_list=[]):
    """Load in the annotation.

    Create a start_dict and tss collection for a set of refseq IDs for a given genome.

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

    mouse_convert_file = os.path.join(annotation_folder, "HMD_HumanPhenotype.rpt")

    # making a dictionary for mouse to human conversion
    mouse_convert_dict = defaultdict(str)

    mouse_convert_table = utils.parse_table(mouse_convert_file, "\t")
    for line in mouse_convert_table:
        mouse_convert_dict[line[4]] = line[0]

    annot_file = genome_dict[genome.upper()]

    start_dict = utils.make_start_dict(annot_file, gene_list)
    tss_loci = []
    if not gene_list:
        gene_list = [*start_dict]
    for gene in gene_list:
        tss_loci.append(utils.make_tss_locus(gene, start_dict, tss_window, tss_window))

    tss_collection = utils.LocusCollection(tss_loci, 50)

    return start_dict, tss_collection, mouse_convert_dict


def split_regions(input_gff, tss_collection, mask_file=None):
    """Split regions if even a single coordinate is shared with the +/-1kb."""
    # create mask regions collection
    if mask_file:
        print("USING MASK FILE {}".format(mask_file))
        # if it's a bed file
        if mask_file.split(".")[-1].upper() == "BED":
            mask_gff = utils.bed_to_gff(mask_file)
        elif mask_file.split(".")[-1].upper() == "GFF":
            mask_gff = utils.parse_table(mask_file, "\t")
        else:
            print("MASK MUST BE A .gff or .bed FILE")

        mask_collection = utils.gff_to_locus_collection(mask_gff)
        print("LOADING {} MASK REGIONS".format(len(mask_collection)))

    split_gff = []
    for line in input_gff:
        chrom = line[0]
        region_id = line[1]
        line_locus = utils.Locus(line[0], line[3], line[4], ".")

        # mask regions
        if mask_file:
            if mask_collection.get_overlap(line_locus, "both"):
                continue

        overlapping_loci = tss_collection.get_overlap(line_locus)
        if overlapping_loci:  # case where a tss overlap
            # identify the parts of the line locus that are contained
            local_tss_collection = utils.LocusCollection(overlapping_loci, 50)
            overlapping_coords = line_locus.coords()
            for tss_locus in overlapping_loci:
                overlapping_coords += tss_locus.coords()

            overlapping_coords = utils.uniquify(overlapping_coords)
            overlapping_coords.sort()

            # you need to hack and slash add 1 to the last coordinate of the overlapping_coords
            overlapping_coords[-1] += 1

            i = 0
            region_ticker = 1
            while i < (len(overlapping_coords) - 1):
                start = int(overlapping_coords[i])
                stop = int(overlapping_coords[(i + 1)]) - 1
                if (stop - start) < 50:  # this eliminates really tiny regions
                    i += 1
                    continue
                split_locus = utils.Locus(chrom, start + 1, stop, ".")

                if line_locus.overlaps(split_locus):
                    new_id = "{}_{}".format(region_id, region_ticker)
                    tss_status = 0
                    if local_tss_collection.get_overlap(split_locus):
                        tss_status = 1
                    split_gff_line = [
                        chrom,
                        new_id,
                        new_id,
                        start,
                        stop,
                        "",
                        ".",
                        tss_status,
                        new_id,
                    ]

                    split_gff.append(split_gff_line)
                    region_ticker += 1
                i += 1
        else:
            line[7] = 0
            split_gff.append(line)

    return split_gff


def map_bams(bam_file_list, split_gff_path, analysis_name, mapped_folder):
    """Map bams to a GFF."""
    print("MAPPING TO THE FOLLOWING BAMS:")

    for bam_file in bam_file_list:
        print(bam_file)
        bam_file_name = os.path.basename(bam_file)

        # MAPPING TO THE STITCHED GFF
        mapped_out_folder = os.path.join(
            mapped_folder, "{}_{}_MAPPED".format(analysis_name, bam_file_name),
        )
        mapped_out_file = os.path.join(mapped_out_folder, "matrix.txt")
        if utils.check_output(mapped_out_file, 0.2, 0.2):
            print(
                "FOUND {} MAPPING DATA FOR BAM: {}".format(
                    split_gff_path, mapped_out_file
                )
            )
        else:
            cmd = "bamliquidator_batch --sense . -e 200 --match_bamToGFF -r {} -o {} {}".format(
                split_gff_path, mapped_out_folder, bam_file,
            )
            print(cmd)

            os.system(cmd)
            if utils.check_output(mapped_out_file, 0.2, 5):
                print(
                    "SUCCESSFULLY MAPPED TO {} FROM BAM: {}".format(
                        split_gff_path, bam_file_name,
                    )
                )
            else:
                print(
                    "ERROR: FAILED TO MAP {} FROM BAM: {}".format(
                        split_gff_path, bam_file_name,
                    )
                )
                sys.exit()

    print("BAM MAPPING COMPLETED NOW MAPPING DATA TO REGIONS")

    # now we make a signal table
    # set up the table using the first bam
    if len(bam_file_list) > 1:

        # set up the first pass at the table
        signal_table = [
            ["REGION_ID", "locusLine"]
            + [name.split("/")[-1] for name in bam_file_list],
        ]
        bam_file_name = bam_file_list[0].split("/")[-1]
        mapped_table = utils.parse_table(
            os.path.join(
                mapped_folder,
                "{}_{}_MAPPED".format(analysis_name, bam_file_name),
                "matrix.txt",
            ),
            "\t",
        )
        signal_table = mapped_table[1:]

        for bam_file in bam_file_list[1:]:
            bam_file_name = bam_file.split("/")[-1]

            mapped_table = utils.parse_table(
                os.path.join(
                    mapped_folder,
                    "{}_{}_MAPPED".format(analysis_name, bam_file_name),
                    "matrix.txt",
                ),
                "\t",
            )

            for i in range(1, len(mapped_table)):
                map_signal = mapped_table[i][2]
                signal_table[i].append(map_signal)
    else:
        bam_file_name = bam_file_list[0].split("/")[-1]
        signal_table = utils.parse_table(
            os.path.join(
                mapped_folder,
                "{}_{}_MAPPED".format(analysis_name, bam_file_name),
                "matrix.txt",
            ),
            "\t",
        )

    return signal_table


def make_average_table(output_folder, analysis_name, use_background=False):
    """Makes a signal table that is the average background subtracted signal for each region.

    If background is present, will zero out regions before trying to take average.
    i.e. no negative regions allowed.

    """
    signal_table_path = os.path.join(
        output_folder, "{}_signal_table.txt".format(analysis_name),
    )
    signal_table = utils.parse_table(signal_table_path, "\t")

    average_table = [["GENE_ID", "locusLine", "{}_signal".format(analysis_name)]]

    # first the easy case with no background
    if not use_background:
        for line in signal_table[1:]:
            new_line = line[0:2]
            avg_signal = round(numpy.mean([float(x) for x in line[2:]]), 4)
            new_line.append(avg_signal)
            average_table.append(new_line)

    # now the condition w/ background
    else:
        control_table_path = os.path.join(
            output_folder, "{}_control_signal_table.txt".format(analysis_name),
        )
        control_table = utils.parse_table(control_table_path, "\t")

        # checking to make sure the # of backgrounds = number of signal bams
        # otherwise throw an error
        signal_n_col = len(signal_table[0])
        control_n_col = len(control_table[0])

        if signal_n_col != control_n_col:
            print("ERROR: MUST PROVIDE SAME NUMBER OF CONTROL BAMS")
            sys.exit()

        signal_n_rows = len(signal_table)
        control_n_rows = len(control_table)

        if signal_n_rows != control_n_rows:
            print("ERROR: MAPPED FILES ARE NOT THE SAME LENGTH")
            sys.exit()

        for i in range(1, len(signal_table)):
            signal_line = signal_table[i]
            control_line = control_table[i]
            if signal_line[0:2] != control_line[0:2]:
                print("ERROR: REGIONS ON LINE {} DO NOT CORRESPOND".format(str(i)))
                sys.exit()

            new_line = signal_line[0:2]

            signal_values = [float(x) for x in signal_line[2:]]
            control_values = [float(x) for x in control_line[2:]]

            subtracted_values = [
                signal_values[x] - control_values[x] for x in range(len(signal_values))
            ]
            subtracted_values = [
                max(0, x) for x in subtracted_values
            ]  # now make negative numbers 0
            avg_signal = round(numpy.mean(subtracted_values), 4)
            new_line.append(avg_signal)
            average_table.append(new_line)

    return average_table


def make_peak_table(
    param_dict,
    split_gff_path,
    average_table_path,
    start_dict,
    gene_list,
    genome_directory,
    tss_window,
    distal_window,
    tads_path="",
):
    """Makes the final peak table with ebox info."""
    peak_table = [
        [
            "REGION_ID",
            "CHROM",
            "START",
            "STOP",
            "LENGTH",
            "TSS",
            "CPG",
            "CPG_FRACTION",
            "GC_FREQ",
            "SIGNAL",
            "CANON_EBOX_COUNT",
            "NON_CANON_EBOX_COUNT",
            "TOTAL_EBOX_COUNT",
            "OVERLAPPING_GENES",
            "PROXIMAL_GENES",
        ]
    ]

    print("LOADING PEAK REGIONS")
    peak_gff = utils.parse_table(split_gff_path, "\t")

    print("LOADING BINDING DATA")
    signal_table = utils.parse_table(average_table_path, "\t")

    print("LOADING CPGS ISLANDS")
    cpg_bed = utils.parse_table(param_dict["cpg_path"], "\t")
    cpg_loci = []
    for line in cpg_bed:
        cpg_loci.append(utils.Locus(line[0], line[1], line[2], ".", line[-1]))
    cpg_collection = utils.LocusCollection(cpg_loci, 50)

    print("MAKING TSS COLLECTIONS")
    if not gene_list:
        gene_list = [*start_dict]

    tss_prox_loci = []
    tss_distal_loci = []
    for ref_id in gene_list:
        tss_prox_loci.append(
            utils.make_tss_locus(ref_id, start_dict, tss_window, tss_window)
        )
        tss_distal_loci.append(
            utils.make_tss_locus(ref_id, start_dict, distal_window, distal_window,)
        )

    # make a 1kb flanking and 50kb flanking collection
    tss_prox_collection = utils.LocusCollection(tss_prox_loci, 50)
    tss_distal_collection = utils.LocusCollection(tss_distal_loci, 50)

    if tads_path:
        print("LOADING TADS FROM {}".format(tads_path))
        tad_collection = utils.import_bound_region(tads_path, "tad")
        use_tads = True

        # building a tad dict keyed by tad ID w/ genes in that tad provided
        tad_dict = defaultdict(list)
        for tss_locus in tss_prox_loci:
            overlapping_tads = tad_collection.get_overlap(tss_locus, "both")
            for tad_locus in overlapping_tads:
                tad_dict[tad_locus.id].append(tss_locus.id)
    else:
        use_tads = False

    print("CLASSIFYING PEAKS")
    ticker = 0

    no_tad_count = 0
    for i in range(len(peak_gff)):
        if not ticker % 1000:
            print(ticker)
        ticker += 1

        # getting the particulars of the region
        gff_line = peak_gff[i]
        peak_id = gff_line[1]
        chrom = gff_line[0]
        start = int(gff_line[3])
        stop = int(gff_line[4])
        line_locus = utils.Locus(chrom, start, stop, ".", peak_id)

        # getting the mapped signal
        signal_line = signal_table[(i + 1)]
        signal_vector = [float(x) for x in signal_line[2:]]

        # setting up the new line
        new_line = [peak_id, chrom, start, stop, line_locus.len()]

        # get the tss status from the gff itself
        # (we are able to do this nicely from the split gff code earlier)
        new_line.append(gff_line[7])

        # check cpg status
        if cpg_collection.get_overlap(line_locus, "both"):
            new_line.append(1)
        else:
            new_line.append(0)

        # now do fractional cpgoverlap
        overlapping_cpg_loci = cpg_collection.get_overlap(line_locus, "both")
        overlapping_bases = 0
        for locus in overlapping_cpg_loci:
            cpg_start = max(locus.start, line_locus.start)
            cpg_end = min(locus.end, line_locus.end)
            overlapping_bases += cpg_end - cpg_start
        overlap_fraction = float(overlapping_bases) / line_locus.len()

        new_line.append(round(overlap_fraction, 2))

        # now get the seq
        line_seq = utils.fetch_seq(genome_directory, chrom, start, stop, True).upper()
        if not line_seq:
            print("UH OH")
            print(line_seq)
            print(gff_line)
            print(i)
            print(chrom)
            print(start)
            print(stop)
            sys.exit()

        gc_freq = float(line_seq.count("GC") + line_seq.count("CG")) / len(line_seq)
        new_line.append(gc_freq)

        # this is where we add the ChIP-seq signal
        new_line += signal_vector

        ebox_match_list = re.findall("CA..TG", line_seq)
        if not ebox_match_list:
            new_line += [0] * 3
        else:
            total_count = len(ebox_match_list)
            canon_count = ebox_match_list.count("CACGTG")
            other_count = total_count - canon_count
            new_line += [canon_count, other_count, total_count]

        # now find the overlapping and proximal genes
        # here each overlapping gene the tss prox locus overlaps the peak

        if use_tads:
            tad_loci = tad_collection.get_overlap(line_locus, "both")

            tad_id_list = [tad_locus.id for tad_locus in tad_loci]
            tad_genes = []
            for tad_id in tad_id_list:
                tad_genes += tad_dict[tad_id]
            if not tad_genes:
                no_tad_count += 1
        else:
            tad_genes = []

        if tad_genes:
            overlapping_genes = [
                start_dict[locus.id]["name"]
                for locus in tss_prox_collection.get_overlap(line_locus, "both")
                if tad_genes.count(locus.id)
            ]
            proximal_genes = [
                start_dict[locus.id]["name"]
                for locus in tss_distal_collection.get_overlap(line_locus, "both")
                if tad_genes.count(locus.id)
            ]
        else:
            overlapping_genes = [
                start_dict[locus.id]["name"]
                for locus in tss_prox_collection.get_overlap(line_locus, "both")
            ]
            proximal_genes = [
                start_dict[locus.id]["name"]
                for locus in tss_distal_collection.get_overlap(line_locus, "both")
            ]

        overlapping_genes = utils.uniquify(overlapping_genes)
        # here the tss 50kb locus overlaps the peak
        # overlap takes priority over proximal
        proximal_genes = [
            gene for gene in proximal_genes if not overlapping_genes.count(gene)
        ]
        proximal_genes = utils.uniquify(proximal_genes)

        overlapping_string = ",".join(overlapping_genes)
        proximal_string = ",".join(proximal_genes)

        new_line += [overlapping_string, proximal_string]

        peak_table.append(new_line)

    print(
        "Out of {} regions, {} were assigned to at least 1 tad".format(
            str(len(peak_table)), str(no_tad_count),
        )
    )
    return peak_table


def make_gene_table(peak_table, analysis_name):
    """Takes the peak table and makes a gene centric table."""
    gene_dict = dict()

    gene_table = [
        [
            "GENE",
            "{}_TSS_SIGNAL".format(analysis_name),
            "{}_DISTAL_SIGNAL".format(analysis_name),
        ],
    ]

    # now iterate through the table
    for line in peak_table[1:]:
        region_length = int(line[4])

        signal = float(line[9]) * region_length

        # genes where this particular peak overlaps the tss prox window
        # where there are both overlap and proximal meet
        if len(line) == 15:
            overlap_gene_list = [gene for gene in line[-2].split(",") if gene]
            if overlap_gene_list.count("107"):
                print(line)
                sys.exit()
            for overlap_gene in overlap_gene_list:
                if overlap_gene not in gene_dict:
                    gene_dict[overlap_gene] = {"tss": 0.0, "distal": 0.0}
                # there can be a nasty 1 overlap case where the region might overlap by the
                # overlapping gene list, but not be real
                if int(line[5]) == 1:
                    gene_dict[overlap_gene]["tss"] += signal
                # this is the case where the gene site is just outside of the promoter or
                # overlapping the gene locus/body these are rare
                else:
                    gene_dict[overlap_gene]["distal"] += signal

            proximal_gene_list = [gene for gene in line[-1].split(",") if gene]
            for proximal_gene in proximal_gene_list:
                if proximal_gene not in gene_dict:
                    gene_dict[proximal_gene] = {"tss": 0.0, "distal": 0.0}
                if not int(line[5]):
                    gene_dict[proximal_gene]["distal"] += signal
        # where there's just overlap
        elif len(line) == 14:
            overlap_gene_list = [gene for gene in line[-1].split(",") if gene]
            if overlap_gene_list.count("107"):
                print(line)
                sys.exit()
            for overlap_gene in overlap_gene_list:
                if overlap_gene not in gene_dict:
                    gene_dict[overlap_gene] = {"tss": 0.0, "distal": 0.0}
                # there can be a nasty 1 overlap case where the region might overlap by the
                # overlapping gene list, but not be real
                if int(line[5]) == 1:
                    gene_dict[overlap_gene]["tss"] += signal
                # this is the case where the mycn site is just outside of the promoter or
                # overlapping the gene locus/body these are rare
                else:
                    gene_dict[overlap_gene]["distal"] += signal

    gene_list = [*gene_dict]
    gene_list = utils.uniquify(gene_list)
    gene_list.sort()

    for gene in gene_list:
        new_line = [gene]
        new_line.append(gene_dict[gene]["tss"])
        new_line.append(gene_dict[gene]["distal"])
        gene_table.append(new_line)

    return gene_table


def call_r_waterfall(gene_table_path, output_folder, analysis_name, top):
    """Function to call the Rscript.

    Wait until the .cls and .gct files are created.
    Returns the paths.

    """
    r_bash_file_path = os.path.join(
        output_folder, "{}_R_plotting.sh".format(analysis_name)
    )
    with open(r_bash_file_path, "w") as r_bash_file:
        r_bash_file.write("#!/usr/bin/bash\n\n")

        r_script_path = os.path.join(
            ROOT_DIR, "scripts", "enhancerPromoter_waterfall.R"
        )
        r_cmd = "Rscript {} {} {} {} {}".format(
            r_script_path,
            gene_table_path,
            "{}/".format(output_folder),
            analysis_name,
            top,
        )
        r_bash_file.write(r_cmd)

    print("writing R plotting command to disk and calling %{}".format(r_bash_file_path))
    os.system("bash {}".format(r_bash_file_path))

    # now check for the .cls output
    cls_path = os.path.join(
        output_folder, "{}_top_{}.cls".format(analysis_name, str(top))
    )

    if utils.check_output(cls_path, 0.5, 5):
        return
    else:
        print(
            "ERROR: UNABLE TO SUCCESFULLY DETECT R SCRIPT OUTPUT AT {}".format(cls_path)
        )
        sys.exit()


# ================================================================================
# ===============================MAIN RUN=========================================
# ================================================================================


def main():
    """Main run method for enhancer promoter contribution tool."""
    parser = argparse.ArgumentParser()

    # required flags
    parser.add_argument(
        "-b",
        "--bam",
        dest="bam",
        nargs="*",
        help="Enter a space separated list of .bam files for the main factor",
        required=True,
    )
    parser.add_argument(
        "-i",
        "--input",
        dest="input",
        type=str,
        help="Enter .gff or .bed file of regions to analyze",
        required=True,
    )
    parser.add_argument(
        "-g",
        "--genome",
        dest="genome",
        type=str,
        help=(
            "specify a genome, HG18,HG19,HG38,MM8,MM9,MM10,RN6 are currently "
            "supported"
        ),
        required=True,
    )
    parser.add_argument(
        "-p",
        "--chrom-path",
        dest="chrom_path",
        type=str,
        help=(
            "Provide path to a folder with a seperate fasta file for each " "chromosome"
        ),
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

    # additional options flags and optional arguments
    parser.add_argument(
        "-a",
        "--activity",
        dest="activity",
        type=str,
        help=(
            "specify a table where first column represents a list of active "
            "refseq genes"
        ),
        required=False,
    )
    parser.add_argument(
        "-c",
        "--control",
        dest="control",
        nargs="*",
        help=(
            "Enter a space separated list of .bam files for background. If "
            "flagged, will perform background subtraction"
        ),
        required=False,
    )
    parser.add_argument(
        "-t",
        "--tss",
        dest="tss",
        type=int,
        help="Define the TSS area +/- the TSS. Default is 1kb",
        required=False,
        default=1000,
    )
    parser.add_argument(
        "-d",
        "--distal",
        dest="distal",
        type=int,
        help="Enter a window to assign distal enhancer signal. Default is 50kb",
        required=False,
        default=50000,
    )
    parser.add_argument(
        "--other-bams",
        dest="other",
        nargs="*",
        help="enter a space separated list of other bams to map to",
        required=False,
    )
    parser.add_argument(
        "--name",
        dest="name",
        type=str,
        help=(
            "enter a root name for the analysis, otherwise will try to find the "
            "name from the input file"
        ),
        required=False,
    )
    parser.add_argument(
        "--top",
        dest="top",
        type=int,
        help=("Run the analysis on the top N genes by total signal. Default is 5000"),
        required=False,
        default=5000,
    )
    parser.add_argument(
        "--tads",
        dest="tads",
        type=str,
        help=("Include a .bed of tad regions to restrict enhancer/gene association"),
        required=False,
        default=None,
    )
    parser.add_argument(
        "--mask",
        dest="mask",
        default=None,
        help=(
            "Mask a set of regions from analysis.  Provide a .bed or .gff of "
            "masking regions"
        ),
    )

    args = parser.parse_args()

    print(args)

    # =====================================================================================
    # ===============================I. PARSING ARGUMENTS==================================
    # =====================================================================================

    print(
        "\n\n#======================================\n#===========I. DATA SUMMARY============\n#="
        "=====================================\n"
    )

    # top analysis subset
    top = args.top

    # input genome
    genome = args.genome.upper()
    print("PERFORMING ANALYSIS ON {} GENOME BUILD".format(genome))

    # set of bams
    bam_file_list = args.bam

    # bring in the input path
    input_path = args.input

    # try to get the input name or use the name argument
    if args.name:
        analysis_name = args.name
    else:
        analysis_name = os.path.basename(input_path).split(".")[0]

    print("USING {} AS ANALYSIS NAME".format(analysis_name))
    # setting up the output folder
    parent_folder = utils.format_folder(args.output, True)
    output_folder = utils.format_folder(
        os.path.join(parent_folder, analysis_name), True
    )

    print("WRITING OUTPUT TO {}".format(output_folder))

    if input_path.split(".")[-1] == "bed":
        # type is bed
        print("input in bed format, converting to gff")
        input_gff = utils.bed_to_gff(input_path)
    else:
        input_gff = utils.parse_table(input_path, "\t")

    # the tss window for proximal signal assignment
    tss_window = int(args.tss)

    # the distal window for assigning nearby enhancer signal
    distal_window = int(args.distal)

    # activity path
    if args.activity:
        activity_path = args.activity
        activity_table = utils.parse_table(activity_path, "\t")
        ref_col = 0
        # try to find the column for refseq id
        for i in range(len(activity_table[2])):  # use an internal row in case of header
            if str(activity_table[1][i]).count("NM_") or str(
                activity_table[1][i]
            ).count("NR_"):
                ref_col = i

        # now check for header
        if not str(activity_table[0][i]).count("NM_") and not str(
            activity_table[0][i]
        ).count("NR_"):
            print("REMOVING HEADER FROM GENE TABLE:")
            print(activity_table[0])
            activity_table.pop(0)

        gene_list = [
            line[ref_col] for line in activity_table
        ]  # this needs to be REFSEQ NM ID
        print("IDENTIFIED {} ACTIVE GENES".format(len(gene_list)))

    else:
        gene_list = []

    # check if tads are being invoked
    if args.tads:
        print("LOADING TAD LOCATIONS FROM {}".format(args.tads))
        tads_path = args.tads
    else:
        tads_path = ""

    print("LOADING ANNOTATION DATA FOR GENOME {}".format(genome))

    genome_dir = args.chrom_path

    # making a chrom_dict that is a list of all chroms with sequence
    chrom_list = utils.uniquify(
        [name.split(".")[0] for name in os.listdir(genome_dir) if name]
    )

    # important here to define the window
    start_dict, tss_collection, mouse_convert_dict = load_annot_file(
        genome, tss_window, gene_list,
    )

    print("FILTERING THE INPUT GFF FOR GOOD CHROMOSOMES")

    print(chrom_list)
    filtered_gff = [line for line in input_gff if chrom_list.count(line[0])]

    print(
        "{} of INITIAL {} REGIONS ARE IN GOOD CHROMOSOMES".format(
            str(len(filtered_gff)), str(len(input_gff)),
        )
    )

    # =====================================================================================
    # ================II. IDENTIFYING TSS PROXIMAL AND DISTAL ELEMENTS=====================
    # =====================================================================================

    print(
        "\n\n#======================================\n#==II. MAPPING TO TSS/DISTAL REGIONS===\n#="
        "=====================================\n"
    )

    # now we need to split the input region
    print("SPLITTING THE INPUT GFF USING A WINDOW OF {}".format(tss_window))
    split_gff = split_regions(filtered_gff, tss_collection, mask_file=args.mask)
    print(len(filtered_gff))
    print(len(split_gff))

    split_gff_path = os.path.join(output_folder, "{}_SPLIT.gff".format(analysis_name))
    utils.unparse_table(split_gff, split_gff_path, "\t")
    print("WRITING TSS SPLIT GFF OUT TO {}".format(split_gff_path))

    # now you have to map the bams to the gff
    print("MAPPING TO THE SPLIT GFF")
    mapped_folder = utils.format_folder(
        os.path.join(output_folder, "bam_mapping"), True
    )

    signal_table = map_bams(bam_file_list, split_gff_path, analysis_name, mapped_folder)
    signal_table_path = os.path.join(
        output_folder, "{}_signal_table.txt".format(analysis_name)
    )
    utils.unparse_table(signal_table, signal_table_path, "\t")

    if args.control:
        control_bam_file_list = args.control
        control_signal_table = map_bams(
            control_bam_file_list, split_gff_path, analysis_name, mapped_folder,
        )
        control_signal_table_path = os.path.join(
            output_folder, "{}_control_signal_table.txt".format(analysis_name),
        )
        utils.unparse_table(control_signal_table, control_signal_table_path, "\t")

    # now create the background subtracted summarized average table
    print("CREATING AN AVERAGE SIGNAL TABLE")
    average_table = make_average_table(
        output_folder, analysis_name, use_background=args.control  # TODO: fix to True or False
    )
    average_table_path = os.path.join(
        output_folder, "{}_average_table.txt".format(analysis_name)
    )
    utils.unparse_table(average_table, average_table_path, "\t")

    # now load up all of the cpg and other parameters to make the actual peak table

    # first check if this has already been done
    peak_table_path = os.path.join(
        output_folder, "{}_PEAK_TABLE.txt".format(analysis_name)
    )
    if utils.check_output(peak_table_path, 0.1, 0.1):
        print("PEAK TABLE OUTPUT ALREADY EXISTS")
        peak_table = utils.parse_table(peak_table_path, "\t")
    else:
        peak_table = make_peak_table(
            param_dict,
            split_gff_path,
            average_table_path,
            start_dict,
            gene_list,
            genome_dir,
            tss_window,
            distal_window,
            tads_path,
        )
        utils.unparse_table(peak_table, peak_table_path, "\t")

    gene_table = make_gene_table(peak_table, analysis_name)

    gene_table_path = os.path.join(
        output_folder, "{}_GENE_TABLE.txt".format(analysis_name)
    )
    utils.unparse_table(gene_table, gene_table_path, "\t")

    # if mouse, need to convert genes over
    if genome.count("MM") == 1:
        print("CONVERTING MOUSE NAMES TO HUMAN HOMOLOGS FOR GSEA")
        converted_gene_table_path = os.path.join(
            output_folder, "{}_GENE_TABLE_CONVERTED.txt".format(analysis_name),
        )

        converted_gene_table = [gene_table[0]]
        for line in gene_table[1:]:
            converted_name = mouse_convert_dict[line[0]]
            if converted_name:
                converted_gene_table.append([converted_name] + line[1:])

                utils.unparse_table(
                    converted_gene_table, converted_gene_table_path, "\t"
                )

        gene_table_path = converted_gene_table_path
        gene_table = converted_gene_table

    # =====================================================================================
    # ===================================III. PLOTTING ====================================
    # =====================================================================================

    print(
        "\n\n#======================================\n#===III. PLOTTING ENHANCER/PROMOTER===\n#=="
        "====================================\n"
    )

    # if there are fewer genes in the gene table than the top genes, only run on all
    if len(gene_table) < int(top):
        print(
            "WARNING: ONLY {} GENES WITH SIGNAL AT EITHER PROMOTERS OR ENHANCERS. NOT ENOUGH TO"
            "RUN ANALYSIS ON TOP {}".format(str(len(gene_table) - 1), str(top))
        )
        top = 0

    # now call the R code
    print("CALLING R PLOTTING SCRIPTS")
    call_r_waterfall(gene_table_path, output_folder, analysis_name, top)


main()
