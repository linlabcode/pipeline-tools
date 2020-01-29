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


# 140306_singleGeneMapper.py
# main method wrapped script to take the enhancer region table output of ROSE_Main
# and map genes to it
# will create two outputs a gene mapped region table where each row is an enhancer
# and a gene table where each row is a gene
# does this by default for super-enhancers only
# update to the gene mapper that finds nearest gene w/ highest signal
# also switching to using the pipeline utils module as opposed to the
# stripped down ROSE_utils module
import argparse
import os
import subprocess
import sys
from collections import defaultdict

from pipeline_tools.definitions import ROOT_DIR
from pipeline_tools.utils import utils

# ==================================================================
# ===========MAPPING GENES TO ENHANCERS WITHOUT BAM RANKING=========
# ==================================================================

# this is the traditional way of running gene mapper


def map_enhancer_to_gene(
    annot_file,
    enhancer_file,
    transcribed_file="",
    unique_genes=True,
    search_window=50000,
    no_format_table=False,
):
    """Maps genes to enhancers.

    If unique_genes, reduces to gene name only. Otherwise, gives for each refseq.

    """
    start_dict = utils.make_start_dict(annot_file)
    enhancer_table = utils.parse_table(enhancer_file, "\t")

    if len(transcribed_file) > 0:
        transcribed_table = utils.parse_table(transcribed_file, "\t")
        # figure out which column has refseq identifiers
        line = transcribed_table[0]
        ref_col = 0
        if len(line) > 1:
            for i in range(len(line)):

                if line[i].count("NM_") > 0 or line[i].count("NR_") > 0:
                    ref_col = i
        transcribed_genes = [line[ref_col] for line in transcribed_table]
    else:
        transcribed_genes = list(start_dict.keys())

    print("MAKING TRANSCRIPT COLLECTION")
    transcribed_collection = utils.make_transcript_collection(
        annot_file, 0, 0, 500, transcribed_genes
    )

    print("MAKING TSS COLLECTION")
    tss_loci = []
    for gene_id in transcribed_genes:
        tss_loci.append(utils.make_tss_locus(gene_id, start_dict, 0, 0))

    # this turns the tss_loci list into a LocusCollection
    # 50 is the internal parameter for LocusCollection and doesn't really matter
    tss_collection = utils.LocusCollection(tss_loci, 50)

    gene_dict = {"overlapping": defaultdict(list), "proximal": defaultdict(list)}

    # dictionaries to hold ranks and super_status of gene nearby enhancers
    rank_dict = defaultdict(list)
    super_dict = defaultdict(list)

    # list of all genes that appear in this analysis
    overall_gene_list = []

    # find the header
    for line in enhancer_table:
        if line[0][0] != "#":
            header = line
            print("this is the header")
            print(header)
            break

    if no_format_table:
        # set up the output tables
        # first by enhancer
        enhancer_to_gene_table = [
            header + ["OVERLAP_GENES", "PROXIMAL_GENES", "CLOSEST_GENE"]
        ]

    else:
        # set up the output tables
        # first by enhancer
        enhancer_to_gene_table = [
            header[0:9]
            + ["OVERLAP_GENES", "PROXIMAL_GENES", "CLOSEST_GENE"]
            + header[-2:]
        ]

    # next make the gene to enhancer table
    gene_to_enhancer_table = [
        [
            "GENE_NAME",
            "REFSEQ_ID",
            "PROXIMAL_ENHANCERS",
            "ENHANCER_RANKS",
            "IS_SUPER",
            "ENHANCER_SIGNAL",
        ]
    ]

    for line in enhancer_table:
        if line[0][0] == "#" or line[0][0] == "R":
            continue

        enhancer_string = "{}:{}-{}".format(line[1], line[2], line[3])

        enhancer_locus = utils.Locus(line[1], line[2], line[3], ".", line[0])

        # overlapping genes are transcribed genes whose transcript is directly in the
        # stitched_locus
        overlapping_loci = transcribed_collection.get_overlap(enhancer_locus, "both")
        overlapping_genes = []
        for overlap_locus in overlapping_loci:
            overlapping_genes.append(overlap_locus.id)

        # proximal_genes are transcribed genes where the tss is within 50kb of the boundary of the
        # stitched loci
        proximal_loci = tss_collection.get_overlap(
            utils.make_search_locus(enhancer_locus, search_window, search_window),
            "both",
        )
        proximal_genes = []
        for prox_locus in proximal_loci:
            proximal_genes.append(prox_locus.id)

        distal_loci = tss_collection.get_overlap(
            utils.make_search_locus(enhancer_locus, 1000000, 1000000), "both"
        )
        distal_genes = []
        for prox_locus in distal_loci:
            distal_genes.append(prox_locus.id)

        overlapping_genes = utils.uniquify(overlapping_genes)
        proximal_genes = utils.uniquify(proximal_genes)
        distal_genes = utils.uniquify(distal_genes)
        all_enhancer_genes = overlapping_genes + proximal_genes + distal_genes
        # these checks make sure each gene list is unique.
        # technically it is possible for a gene to be overlapping, but not proximal since the
        # gene could be longer than the 50kb window, but we'll let that slide here
        for ref_id in overlapping_genes:
            if proximal_genes.count(ref_id) == 1:
                proximal_genes.remove(ref_id)

        for ref_id in proximal_genes:
            if distal_genes.count(ref_id) == 1:
                distal_genes.remove(ref_id)

        # Now find the closest gene
        if len(all_enhancer_genes) == 0:
            closest_gene = ""
        else:
            # get enhancer_center
            enhancer_center = (int(line[2]) + int(line[3])) / 2

            # get absolute distance to enhancer center
            dist_list = [
                abs(enhancer_center - start_dict[gene_id]["start"][0])
                for gene_id in all_enhancer_genes
            ]
            # get the ID and convert to name
            closest_gene = start_dict[
                all_enhancer_genes[dist_list.index(min(dist_list))]
            ]["name"]

        # NOW WRITE THE ROW FOR THE ENHANCER TABLE
        if no_format_table:

            new_enhancer_line = list(line)
            new_enhancer_line.append(
                ",".join(
                    utils.uniquify([start_dict[x]["name"] for x in overlapping_genes])
                )
            )
            new_enhancer_line.append(
                ".".join(
                    utils.uniquify([start_dict[x]["name"] for x in proximal_genes])
                )
            )
            new_enhancer_line.append(closest_gene)

        else:
            new_enhancer_line = line[0:9]
            new_enhancer_line.append(
                ",".join(
                    utils.uniquify([start_dict[x]["name"] for x in overlapping_genes])
                )
            )
            new_enhancer_line.append(
                ",".join(
                    utils.uniquify([start_dict[x]["name"] for x in proximal_genes])
                )
            )
            new_enhancer_line.append(closest_gene)
            new_enhancer_line += line[-2:]

        enhancer_to_gene_table.append(new_enhancer_line)
        # Now grab all overlapping and proximal genes for the gene ordered table

        overall_gene_list += overlapping_genes
        for ref_id in overlapping_genes:
            gene_dict["overlapping"][ref_id].append(enhancer_string)
            rank_dict[ref_id].append(int(line[-2]))
            super_dict[ref_id].append(int(line[-1]))

        overall_gene_list += proximal_genes
        for ref_id in proximal_genes:
            gene_dict["proximal"][ref_id].append(enhancer_string)
            rank_dict[ref_id].append(int(line[-2]))
            super_dict[ref_id].append(int(line[-1]))

    # End loop through

    # Make table by gene
    overall_gene_list = utils.uniquify(overall_gene_list)

    # use enhancer rank to order
    rank_order = utils.order([min(rank_dict[x]) for x in overall_gene_list])

    used_names = []
    for i in rank_order:
        ref_id = overall_gene_list[i]
        gene_name = start_dict[ref_id]["name"]
        if used_names.count(gene_name) > 0 and unique_genes is True:
            continue
        else:
            used_names.append(gene_name)

        prox_enhancers = (
            gene_dict["overlapping"][ref_id] + gene_dict["proximal"][ref_id]
        )

        super_status = max(super_dict[ref_id])
        enhancer_ranks = ",".join([str(x) for x in rank_dict[ref_id]])

        new_line = [
            gene_name,
            ref_id,
            ",".join(prox_enhancers),
            enhancer_ranks,
            super_status,
        ]
        gene_to_enhancer_table.append(new_line)

    # resort enhancer_to_gene_table
    if no_format_table:
        return enhancer_to_gene_table, gene_to_enhancer_table
    else:
        enhancer_order = utils.order(
            [int(line[-2]) for line in enhancer_to_gene_table[1:]]
        )
        sorted_table = [enhancer_to_gene_table[0]]
        for i in enhancer_order:
            sorted_table.append(enhancer_to_gene_table[(i + 1)])

        return sorted_table, gene_to_enhancer_table


# ==================================================================
# ===========MAPPING GENES TO ENHANCERS WITH BAM RANKING============
# ==================================================================


def make_signal_dict(mapped_gff_file, control_mapped_gff_file=""):
    """
    makes a signal dict
    """
    print(
        (
            "\t called make_signal_dict on {} (ctrl: {})".format(
                mapped_gff_file, control_mapped_gff_file
            )
        )
    )
    signal_dict = defaultdict(float)

    mapped_gff = utils.parse_table(mapped_gff_file, "\t")
    if len(control_mapped_gff_file) > 0:
        control_gff = utils.parse_table(control_mapped_gff_file, "\t")

        for i in range(1, len(mapped_gff)):

            signal = float(mapped_gff[i][2]) - float(control_gff[i][2])
            if signal < 0:
                signal = 0.0
            signal_dict[mapped_gff[i][0]] = signal
    else:
        for i in range(1, len(mapped_gff)):
            signal = float(mapped_gff[i][2])
            signal_dict[mapped_gff[i][0]] = signal

    return signal_dict


def map_enhancer_to_gene_top(
    rank_by_bam_file,
    control_bam_file,
    genome,
    annot_file,
    enhancer_file,
    transcribed_file="",
    unique_genes=True,
    search_window=50000,
    no_format_table=False,
):
    """Maps genes to enhancers.

    If unique_genes, reduces to gene name only. Otherwise, gives for each refseq.

    """
    start_dict = utils.make_start_dict(annot_file)
    enhancer_name = enhancer_file.split("/")[-1].split(".")[0]
    enhancer_table = utils.parse_table(enhancer_file, "\t")

    if len(transcribed_file) > 0:
        transcribed_table = utils.parse_table(transcribed_file, "\t")
        transcribed_genes = [line[1] for line in transcribed_table]
    else:
        transcribed_genes = list(start_dict.keys())

    print("MAKING TRANSCRIPT COLLECTION")
    transcribed_collection = utils.make_transcript_collection(
        annot_file, 0, 0, 500, transcribed_genes
    )

    print("MAKING TSS COLLECTION")
    tss_loci = []
    for gene_id in transcribed_genes:
        tss_loci.append(utils.make_tss_locus(gene_id, start_dict, 0, 0))

    # this turns the tss_loci list into a LocusCollection
    # 50 is the internal parameter for LocusCollection and doesn't really
    # matter
    tss_collection = utils.LocusCollection(tss_loci, 50)

    gene_dict = {"overlapping": defaultdict(list), "proximal": defaultdict(list)}

    # dictionaries to hold ranks and super_status of gene nearby enhancers
    rank_dict = defaultdict(list)
    super_dict = defaultdict(list)

    # list of all genes that appear in this analysis
    overall_gene_list = []

    # find the header
    for line in enhancer_table:
        if line[0][0] != "#":
            header = line
            print("this is the header")
            print(header)
            break

    if no_format_table:
        # set up the output tables
        # first by enhancer
        enhancer_to_gene_table = [
            header + ["OVERLAP_GENES", "PROXIMAL_GENES", "CLOSEST_GENE"]
        ]

    else:
        # set up the output tables
        # first by enhancer
        enhancer_to_gene_table = [
            header[0:9]
            + ["OVERLAP_GENES", "PROXIMAL_GENES", "CLOSEST_GENE"]
            + header[-2:]
        ]

    # next make the gene to enhancer table
    gene_to_enhancer_table = [
        [
            "GENE_NAME",
            "REFSEQ_ID",
            "PROXIMAL_ENHANCERS",
            "ENHANCER_RANKS",
            "IS_SUPER",
            "TSS_SIGNAL",
        ]
    ]

    for line in enhancer_table:
        if line[0][0] == "#" or line[0][0] == "R":
            continue

        enhancer_string = "{}:{}-{}".format(line[1], line[2], line[3])

        enhancer_locus = utils.Locus(line[1], line[2], line[3], ".", line[0])

        # overlapping genes are transcribed genes whose transcript is directly
        # in the stitchedLocus
        overlapping_loci = transcribed_collection.get_overlap(enhancer_locus, "both")
        overlapping_genes = []
        for overlap_locus in overlapping_loci:
            overlapping_genes.append(overlap_locus.id)

        # proximal_genes are transcribed genes where the tss is within 50kb of
        # the boundary of the stitched loci
        proximal_loci = tss_collection.get_overlap(
            utils.make_search_locus(enhancer_locus, search_window, search_window),
            "both",
        )
        proximal_genes = []
        for prox_locus in proximal_loci:
            proximal_genes.append(prox_locus.id)

        distal_loci = tss_collection.get_overlap(
            utils.make_search_locus(enhancer_locus, 1000000, 1000000), "both"
        )
        distal_genes = []
        for prox_locus in distal_loci:
            distal_genes.append(prox_locus.id)

        overlapping_genes = utils.uniquify(overlapping_genes)
        proximal_genes = utils.uniquify(proximal_genes)
        distal_genes = utils.uniquify(distal_genes)
        all_enhancer_genes = overlapping_genes + proximal_genes + distal_genes
        # these checks make sure each gene list is unique.
        # technically it is possible for a gene to be overlapping, but not proximal since the
        # gene could be longer than the 50kb window, but we'll let that slide
        # here
        for ref_id in overlapping_genes:
            if proximal_genes.count(ref_id) == 1:
                proximal_genes.remove(ref_id)

        for ref_id in proximal_genes:
            if distal_genes.count(ref_id) == 1:
                distal_genes.remove(ref_id)

        # Now find the closest gene
        if len(all_enhancer_genes) == 0:
            closest_gene = ""
        else:
            # get enhancer_center
            enhancer_center = (int(line[2]) + int(line[3])) // 2

            # get absolute distance to enhancer center
            dist_list = [
                abs(enhancer_center - start_dict[gene_id]["start"][0])
                for gene_id in all_enhancer_genes
            ]
            # get the ID and convert to name
            closest_gene = start_dict[
                all_enhancer_genes[dist_list.index(min(dist_list))]
            ]["name"]

        # NOW WRITE THE ROW FOR THE ENHANCER TABLE
        if no_format_table:

            new_enhancer_line = list(line)
            new_enhancer_line.append(
                ",".join(
                    utils.uniquify([start_dict[x]["name"] for x in overlapping_genes])
                )
            )
            new_enhancer_line.append(
                ",".join(
                    utils.uniquify([start_dict[x]["name"] for x in proximal_genes])
                )
            )
            new_enhancer_line.append(closest_gene)

        else:
            new_enhancer_line = line[0:9]
            new_enhancer_line.append(
                ",".join(
                    utils.uniquify([start_dict[x]["name"] for x in overlapping_genes])
                )
            )
            new_enhancer_line.append(
                ",".join(
                    utils.uniquify([start_dict[x]["name"] for x in proximal_genes])
                )
            )
            new_enhancer_line.append(closest_gene)
            new_enhancer_line += line[-2:]

        enhancer_to_gene_table.append(new_enhancer_line)
        # Now grab all overlapping and proximal genes for the gene ordered
        # table

        overall_gene_list += overlapping_genes
        for ref_id in overlapping_genes:
            gene_dict["overlapping"][ref_id].append(enhancer_string)
            rank_dict[ref_id].append(int(line[-2]))
            super_dict[ref_id].append(int(line[-1]))

        overall_gene_list += proximal_genes
        for ref_id in proximal_genes:
            gene_dict["proximal"][ref_id].append(enhancer_string)
            rank_dict[ref_id].append(int(line[-2]))
            super_dict[ref_id].append(int(line[-1]))

    # End loop through
    # Make table by gene
    print("MAKING ENHANCER ASSOCIATED GENE TSS COLLECTION")
    overall_gene_list = utils.uniquify(overall_gene_list)

    # get the chromLists from the various bams here
    cmd = "samtools idxstats {}".format(rank_by_bam_file)
    idx_stats = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)
    idx_stats = idx_stats.communicate()
    bam_chrom_list = [
        line.split("\t")[0] for line in idx_stats[0].decode("utf-8").split("\n")[0:-2]
    ]

    if len(control_bam_file) > 0:
        cmd = "samtools idxstats {}".format(control_bam_file)
        idx_stats = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)
        idx_stats = idx_stats.communicate()
        bam_chrom_list_control = [
            line.split("\t")[0]
            for line in idx_stats[0].decode("utf-8").split("\n")[0:-2]
        ]
        bam_chrom_list = [
            chrom
            for chrom in bam_chrom_list
            if bam_chrom_list_control.count(chrom) != 0
        ]

    # now make sure no genes have a bad chrom
    overall_gene_list = [
        gene
        for gene in overall_gene_list
        if bam_chrom_list.count(start_dict[gene]["chr"]) != 0
    ]

    # now make an enhancer collection of all transcripts
    enhancer_gene_collection = utils.make_transcript_collection(
        annot_file, 5000, 5000, 500, overall_gene_list
    )

    enhancer_gene_gff = utils.locus_collection_to_gff(enhancer_gene_collection)

    # dump the gff to file
    enhancer_folder = os.path.dirname(enhancer_file)
    gff_root_name = "{}_TSS_ENHANCER_GENES_-5000_+5000".format(genome)
    enhancer_gene_gff_file = os.path.join(
        enhancer_folder, "{}_{}.gff".format(enhancer_name, gff_root_name)
    )
    utils.unparse_table(enhancer_gene_gff, enhancer_gene_gff_file, "\t")

    # now we need to run bam_to_gff

    print("MAPPING SIGNAL AT ENHANCER ASSOCIATED GENE TSS")
    # map density at genes in the +/- 5kb tss region
    # first on the rankBy bam
    bam_name = os.path.basename(rank_by_bam_file)
    mapped_rank_by_folder = os.path.join(
        enhancer_folder, "{}_{}_{}".format(enhancer_name, gff_root_name, bam_name)
    )
    mapped_rank_by_file = os.path.join(
        enhancer_folder,
        "{}_{}_{}".format(enhancer_name, gff_root_name, bam_name),
        "matrix.txt",
    )
    cmd = "bamliquidator_batch --sense . -e 200 --match_bamToGFF -r {} -o {} {}".format(
        enhancer_gene_gff_file, mapped_rank_by_folder, rank_by_bam_file
    )
    print("Mapping rankby bam {}".format(rank_by_bam_file))
    print(cmd)
    os.system(cmd)

    # check for completion
    if utils.check_output(mapped_rank_by_file, 0.2, 5):
        print(
            "SUCCESSFULLY MAPPED TO {} FROM BAM: {}".format(
                enhancer_gene_gff_file, rank_by_bam_file
            )
        )
    else:
        print(
            "ERROR: FAILED TO MAP {} FROM BAM: {}".format(
                enhancer_gene_gff_file, rank_by_bam_file
            )
        )
        sys.exit()

    # next on the control bam if it exists
    if len(control_bam_file) > 0:
        control_name = os.path.basename(control_bam_file)
        mapped_control_folder = os.path.join(
            enhancer_folder,
            "{}_{}_{}".format(enhancer_name, gff_root_name, control_name),
        )
        mapped_control_file = os.path.join(
            enhancer_folder,
            "{}_{}_{}".format(enhancer_name, gff_root_name, control_name),
            "matrix.txt",
        )
        cmd = "bamliquidator_batch --sense . -e 200 --match_bamToGFF -r {} -o {} {}".format(
            enhancer_gene_gff_file, mapped_control_folder, control_bam_file
        )
        print("Mapping control bam {}".format(control_bam_file))
        print(cmd)
        os.system(cmd)

        # check for completion
        if utils.check_output(mapped_control_file, 0.2, 5):
            print(
                "SUCCESSFULLY MAPPED TO {} FROM BAM: {}".format(
                    enhancer_gene_gff_file, control_bam_file
                )
            )
        else:
            print(
                "ERROR: FAILED TO MAP {} FROM BAM: {}".format(
                    enhancer_gene_gff_file, control_bam_file
                )
            )
            sys.exit()

    # now get the appropriate output files
    if len(control_bam_file) > 0:
        print(
            "CHECKING FOR MAPPED OUTPUT AT {} AND {}".format(
                mapped_rank_by_file, mapped_control_file
            )
        )
        if utils.check_output(mapped_rank_by_file, 1, 1) and utils.check_output(
            mapped_control_file, 1, 1
        ):
            print("MAKING ENHANCER ASSOCIATED GENE TSS SIGNAL DICTIONARIES")
            signal_dict = make_signal_dict(mapped_rank_by_file, mapped_control_file)
        else:
            print("NO MAPPING OUTPUT DETECTED")
            sys.exit()
    else:
        print("CHECKING FOR MAPPED OUTPUT AT {}".format(mapped_rank_by_file))
        if utils.check_output(mapped_rank_by_file, 1, 30):
            print("MAKING ENHANCER ASSOCIATED GENE TSS SIGNAL DICTIONARIES")
            signal_dict = make_signal_dict(mapped_rank_by_file)
        else:
            print("NO MAPPING OUTPUT DETECTED")
            sys.exit()

    # use enhancer rank to order

    rank_order = utils.order([min(rank_dict[x]) for x in overall_gene_list])

    used_names = []

    # make a new dict to hold TSS signal by max per gene_name
    gene_name_sig_dict = defaultdict(list)
    print("MAKING GENE TABLE")
    for i in rank_order:
        ref_id = overall_gene_list[i]
        gene_name = start_dict[ref_id]["name"]
        if used_names.count(gene_name) > 0 and unique_genes is True:
            continue
        else:
            used_names.append(gene_name)

        prox_enhancers = (
            gene_dict["overlapping"][ref_id] + gene_dict["proximal"][ref_id]
        )

        super_status = max(super_dict[ref_id])
        enhancer_ranks = ",".join([str(x) for x in rank_dict[ref_id]])

        enhancer_signal = signal_dict[ref_id]
        gene_name_sig_dict[gene_name].append(enhancer_signal)

        new_line = [
            gene_name,
            ref_id,
            ",".join(prox_enhancers),
            enhancer_ranks,
            super_status,
            enhancer_signal,
        ]
        gene_to_enhancer_table.append(new_line)

    print("MAKING ENHANCER TO TOP GENE TABLE")

    if no_format_table:
        enhancer_to_top_gene_table = [
            enhancer_to_gene_table[0] + ["TOP_GENE", "TSS_SIGNAL"]
        ]
    else:
        enhancer_to_top_gene_table = [
            enhancer_to_gene_table[0][0:12]
            + ["TOP_GENE", "TSS_SIGNAL"]
            + enhancer_to_gene_table[0][-2:]
        ]

    for line in enhancer_to_gene_table[1:]:

        gene_list = []
        if no_format_table:
            gene_list += line[-3].split(",")
            gene_list += line[-2].split(",")

        else:
            gene_list += line[10].split(",")
            gene_list += line[11].split(",")

        gene_list = utils.uniquify([x for x in gene_list if len(x) > 0])
        if len(gene_list) > 0:
            try:
                sig_vector = [max(gene_name_sig_dict[x]) for x in gene_list]
                max_index = sig_vector.index(max(sig_vector))
                max_gene = gene_list[max_index]
                max_sig = sig_vector[max_index]
                if max_sig == 0.0:
                    max_gene = "NONE"
                    max_sig = "NONE"
            except ValueError:
                if len(gene_list) == 1:
                    max_gene = gene_list[0]
                    max_sig = "NONE"
                else:
                    max_gene = "NONE"
                    max_sig = "NONE"
        else:
            max_gene = "NONE"
            max_sig = "NONE"
        if no_format_table:
            new_line = line + [max_gene, max_sig]
        else:
            new_line = line[0:12] + [max_gene, max_sig] + line[-2:]
        enhancer_to_top_gene_table.append(new_line)

    # resort enhancer_to_gene_table
    if no_format_table:
        return (
            enhancer_to_gene_table,
            enhancer_to_top_gene_table,
            gene_to_enhancer_table,
        )
    else:
        enhancer_order = utils.order(
            [int(line[-2]) for line in enhancer_to_gene_table[1:]]
        )
        sorted_table = [enhancer_to_gene_table[0]]
        sorted_top_gene_table = [enhancer_to_top_gene_table[0]]
        for i in enhancer_order:
            sorted_table.append(enhancer_to_gene_table[(i + 1)])
            sorted_top_gene_table.append(enhancer_to_top_gene_table[(i + 1)])

        return sorted_table, sorted_top_gene_table, gene_to_enhancer_table


# ==================================================================
# =========================MAIN METHOD==============================
# ==================================================================
def main():
    """Main run call."""
    parser = argparse.ArgumentParser()
    # required flags
    parser.add_argument(
        "-i",
        "--i",
        dest="input",
        default=None,
        required=True,
        help="Enter a ROSE ranked enhancer or super-enhancer file",
    )
    parser.add_argument(
        "-g",
        "--genome",
        dest="genome",
        default=None,
        required=True,
        help="Enter the genome build (MM9,MM8,HG18,HG19)",
    )

    # optional flags
    parser.add_argument(
        "-r",
        "--rankby",
        dest="rankby",
        required=False,
        help="Enter the bam used to rank enhancers",
    )
    parser.add_argument(
        "-c",
        "--control",
        dest="control",
        default='',
        help="Enter a background bam for background correction",
    )
    parser.add_argument(
        "-l",
        "--list",
        dest="gene_list",
        required=False,
        help="Enter a gene list to filter through",
    )
    parser.add_argument(
        "-o",
        "--out",
        dest="out",
        required=False,
        help="Enter an output folder. Default will be same folder as input file",
    )
    parser.add_argument(
        "-w",
        "--window",
        dest="window",
        default=50000,
        help="Enter a search distance for genes. Default is 50,000bp",
    )
    parser.add_argument(
        "-f",
        "--format",
        dest="format_table",
        action="store_true",
        default=False,
        help="If flagged, maintains original formatting of input table",
    )

    args = parser.parse_args()

    print(args)

    # GETTING THE GENOME
    genome = args.genome
    print("USING {} AS THE GENOME".format(genome))

    # GETTING THE CORRECT ANNOT FILE
    annotation_folder = os.path.join(ROOT_DIR, "annotation")
    genome_dict = {
        "HG18": os.path.join(annotation_folder, "hg18_refseq.ucsc"),
        "MM9": os.path.join(annotation_folder, "mm9_refseq.ucsc"),
        "MM10": os.path.join(annotation_folder, "mm10_refseq.ucsc"),
        "HG19": os.path.join(annotation_folder, "hg19_refseq.ucsc"),
        "RN4": os.path.join(annotation_folder, "rn4_refseq.ucsc"),
        "RN6": os.path.join(annotation_folder, "rn6_refseq.ucsc"),
        "HG38": os.path.join(annotation_folder, "hg38_refseq.ucsc"),
    }
    annot_file = genome_dict[genome.upper()]

    # GETTING THE INPUT
    enhancer_file = args.input
    window = int(args.window)

    # making the out folder if it doesn't exist
    if args.out:
        out_folder = utils.format_folder(args.out, True)
    else:
        out_folder = os.path.dirname(enhancer_file)

    # GETTING BAM INFO
    rank_by_bam_file = args.rankby
    control_bam_file = args.control

    # CHECK FORMATTING FLAG
    if args.format_table:
        no_format_table = True
    else:
        no_format_table = False

    # GETTING THE TRANSCRIBED LIST
    if args.gene_list:

        transcribed_file = args.gene_list
    else:
        transcribed_file = ""

    # Writing enhancer output
    enhancer_file_name = os.path.basename(enhancer_file).split(".")[0]
    if args.rankby:
        print("MAPPING GENES TO ENHANCERS USING CLOSEST ACTIVE GENE")
        (
            enhancer_to_gene_table,
            enhancer_to_top_gene_table,
            gene_to_enhancer_table,
        ) = map_enhancer_to_gene_top(
            rank_by_bam_file,
            control_bam_file,
            genome,
            annot_file,
            enhancer_file,
            transcribed_file,
            True,
            window,
            no_format_table,
        )

        if window != 50000:
            # writing the enhancer table

            out1 = os.path.join(
                out_folder,
                "{}_ENHANCER_TO_GENE_{}KB.txt".format(
                    enhancer_file_name, str(window // 1000)
                ),
            )
            print("writing output to {}".format(out1))
            utils.unparse_table(enhancer_to_gene_table, out1, "\t")

            # writing enhancer top gene table
            out2 = os.path.join(
                out_folder,
                "{}_ENHANCER_TO_TOP_GENE_{}KB.txt".format(
                    enhancer_file_name, str(window // 1000)
                ),
            )
            utils.unparse_table(enhancer_to_top_gene_table, out2, "\t")

            # writing the gene table
            out3 = os.path.join(
                out_folder,
                "{}_GENE_TO_ENHANCER_{}KB.txt".format(
                    enhancer_file_name, str(window // 1000)
                ),
            )
            utils.unparse_table(gene_to_enhancer_table, out3, "\t")
        else:
            # writing the enhancer table
            out1 = os.path.join(
                out_folder, "{}_ENHANCER_TO_GENE.txt".format(enhancer_file_name)
            )
            utils.unparse_table(enhancer_to_gene_table, out1, "\t")

            # writing the enhancer table
            out2 = os.path.join(
                out_folder, "{}_ENHANCER_TO_TOP_GENE.txt".format(enhancer_file_name)
            )
            utils.unparse_table(enhancer_to_top_gene_table, out2, "\t")

            # writing the gene table
            out3 = os.path.join(
                out_folder, "{}_GENE_TO_ENHANCER.txt".format(enhancer_file_name)
            )
            utils.unparse_table(gene_to_enhancer_table, out3, "\t")
    else:
        # do traditional mapping
        print("MAPPING GENES TO ENHANCERS USING PROXIMITY RULE")
        enhancer_to_gene_table, gene_to_enhancer_table = map_enhancer_to_gene(
            annot_file, enhancer_file, transcribed_file, True, window, no_format_table
        )

        if window != 50000:
            # writing the enhancer table
            out1 = os.path.join(
                out_folder,
                "{}_ENHANCER_TO_GENE_{}KB.txt".format(
                    enhancer_file_name, str(window // 1000)
                ),
            )
            utils.unparse_table(enhancer_to_gene_table, out1, "\t")

            # writing the gene table
            out2 = os.path.join(
                out_folder,
                "{}_GENE_TO_ENHANCER_{}KB.txt".format(
                    enhancer_file_name, str(window // 1000)
                ),
            )
            utils.unparse_table(gene_to_enhancer_table, out2, "\t")
        else:
            # writing the enhancer table
            out1 = os.path.join(
                out_folder, "{}_ENHANCER_TO_GENE.txt".format(enhancer_file_name)
            )
            utils.unparse_table(enhancer_to_gene_table, out1, "\t")

            # writing the gene table
            out2 = os.path.join(
                out_folder, "{}_GENE_TO_ENHANCER.txt".format(enhancer_file_name)
            )
            utils.unparse_table(gene_to_enhancer_table, out2, "\t")


if __name__ == "__main__":
    main()
