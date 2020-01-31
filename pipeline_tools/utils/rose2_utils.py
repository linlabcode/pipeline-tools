"""ROSE2 utility functions."""

import copy
import os
import subprocess
import sys
from collections import defaultdict

import numpy

from pipeline_tools.definitions import ROOT_DIR
from pipeline_tools.utils import utils


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


# ==================================================================
# =====================HELPER FUNCTIONS=============================
# ==================================================================


def get_bam_chrom_list(bam_file_list):
    """Get the consensus list of chromosomes mapped by the bams."""
    # start w/ the first bam
    cmd = "samtools idxstats {}".format(bam_file_list[0])
    idx_stats = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)
    idx_stats = idx_stats.communicate()
    final_chrom_list = [
        line.split("\t")[0] for line in idx_stats[0].decode("utf-8").split("\n")[:-2]
    ]

    # now go through each additional bam
    for bam_file in bam_file_list:
        cmd = "samtools idxstats {}".format(bam_file)
        idx_stats = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)
        idx_stats = idx_stats.communicate()
        chrom_list = [
            line.split("\t")[0]
            for line in idx_stats[0].decode("utf-8").split("\n")[:-2]
        ]
        final_chrom_list = [
            chrom for chrom in final_chrom_list if chrom_list.count(chrom) != 0
        ]

    return utils.uniquify(final_chrom_list)


def check_ref_collection(reference_collection):
    """Make sure the names of all loci in the reference collection are unique."""
    names_list = [locus.id for locus in reference_collection.get_loci()]
    if len(names_list) != len(utils.uniquify(names_list)):
        print("ERROR: REGIONS HAVE NON-UNIQUE IDENTIFIERS")
        sys.exit()
    else:
        print("REFERENCE COLLECTION PASSES QC")
        return


def filter_gff(gff_file, chrom_list):
    """Take in a gff and filter out all lines that don't belong to a chrom in the chrom_list."""
    gff = utils.parse_table(gff_file, "\t")
    filtered_gff = []
    exclude_list = []
    for line in gff:
        if chrom_list.count(line[0]) == 1:
            filtered_gff.append(line)
        else:
            exclude_list.append(line[0])

    exclude_list = utils.uniquify(exclude_list)
    if len(exclude_list) > 0:
        print(
            "EXCLUDED GFF REGIONS FROM THE FALLING CHROMS: {}".format(
                ",".join(exclude_list)
            )
        )

    return filtered_gff


# ==================================================================
# =====================REGION STITCHING=============================
# ==================================================================


def optimize_stitching(locus_collection, name, out_folder, step_size=500):
    """
    takes a locus collection and starts writing out stitching stats at step sized intervals
    """
    max_stitch = 15000  # set a hard wired match stitching parameter

    stitch_table = [
        [
            "STEP",
            "NUM_REGIONS",
            "TOTAL_CONSTIT",
            "TOTAL_REGION",
            "MEAN_CONSTIT",
            "MEDIAN_CONSTIT",
            "MEAN_REGION",
            "MEDIAN_REGION",
            "MEAN_STITCH_FRACTION",
            "MEDIAN_STITCH_FRACTION",
        ]
    ]
    # first consolidate the collection
    locus_collection = locus_collection.stitch_collection(stitch_window=0)
    total_constit = sum([locus.len() for locus in locus_collection.get_loci()])
    step = 0
    while step <= max_stitch:

        print("Getting stitch stats for {} (bp)".format(step))
        stitch_collection = locus_collection.stitch_collection(stitch_window=step)
        num_regions = len(stitch_collection)
        stitch_loci = stitch_collection.get_loci()
        region_lengths = [locus.len() for locus in stitch_loci]
        total_region = sum(region_lengths)
        constit_lengths = []
        for locus in stitch_loci:
            constit_loci = locus_collection.get_overlap(locus)
            constit_lengths.append(sum([locus.len() for locus in constit_loci]))

        mean_constit = round(numpy.mean(constit_lengths), 2)
        median_constit = round(numpy.median(constit_lengths), 2)

        mean_region = round(numpy.mean(region_lengths), 2)
        median_region = round(numpy.median(region_lengths), 2)

        stitch_fractions = [
            float(constit_lengths[i]) / float(region_lengths[i])
            for i in range(len(region_lengths))
        ]
        mean_stitch_fraction = round(numpy.mean(stitch_fractions), 2)
        median_stitch_fraction = round(numpy.median(stitch_fractions), 2)

        new_line = [
            step,
            num_regions,
            total_constit,
            total_region,
            mean_constit,
            median_constit,
            mean_region,
            median_region,
            mean_stitch_fraction,
            median_stitch_fraction,
        ]

        stitch_table.append(new_line)

        step += step_size

    # write the stitch table to disk
    stitch_param_file = os.path.join(out_folder, "{}_stitch_params.tmp".format(name))
    utils.unparse_table(stitch_table, stitch_param_file, "\t")
    # call the rscript
    r_cmd = "Rscript {} {} {} {}".format(
        os.path.join(ROOT_DIR, "scripts", "ROSE2_stitchOpt.R"),
        stitch_param_file,
        out_folder + "/",  # TODO: fix R script so it does not require '/'
        name,
    )
    print(r_cmd)
    # get back the stitch parameter
    r_output = subprocess.Popen(r_cmd, stdout=subprocess.PIPE, shell=True)
    r_output_test = r_output.communicate()

    print(r_output_test)

    stitch_param = r_output_test[0].decode("utf-8").split("\n")[2]
    try:
        stitch_param = int(stitch_param)
    except ValueError:
        print("INVALID STITCHING PARAMETER. STITCHING OPTIMIZATION FAILED")
        sys.exit()

    # delete? the table
    # os.system('rm -f %s' % (stitch_param_file))
    return stitch_param


def region_stitching(
    reference_collection,
    name,
    out_folder,
    stitch_window,
    tss_window,
    annot_file,
    remove_tss=True,
):
    """Preform region stitching."""
    print("PERFORMING REGION STITCHING")
    # first have to turn bound region file into a locus collection

    # need to make sure this names correctly... each region should have a unique name
    # reference_collection

    debug_output = []
    # filter out all bound regions that overlap the TSS of an ACTIVE GENE
    if remove_tss:

        print(
            "REMOVING TSS FROM REGIONS USING AN EXCLUSION WINDOW OF {}BP".format(
                str(tss_window)
            )
        )
        # first make a locus collection of TSS

        start_dict = utils.make_start_dict(annot_file)

        # now makeTSS loci for active genes
        remove_ticker = 0
        # this loop makes a locus centered around +/- tss_window of transcribed genes
        # then adds it to the list tss_loci
        tss_loci = []
        for gene_id in list(start_dict.keys()):
            tss_loci.append(
                utils.make_tss_locus(gene_id, start_dict, tss_window, tss_window)
            )

        # this turns the tss_loci list into a LocusCollection
        # 50 is the internal parameter for LocusCollection and doesn't really matter
        tss_collection = utils.LocusCollection(tss_loci, 50)

        # gives all the loci in reference_collection
        bound_loci = list(reference_collection.get_loci())

        # this loop will check if each bound region is contained by the TSS exclusion zone
        # this will drop out a lot of the promoter only regions that are tiny
        # typical exclusion window is around 2kb
        for locus in bound_loci:
            if len(tss_collection.get_containers(locus, "both")) > 0:

                # if true, the bound locus overlaps an active gene
                reference_collection.remove(locus)
                debug_output.append([locus.__str__(), locus.id, "CONTAINED"])
                remove_ticker += 1
        print(
            "REMOVED {} LOCI BECAUSE THEY WERE CONTAINED BY A TSS".format(
                str(remove_ticker)
            )
        )

    # reference_collection is now all enriched region loci that don't overlap an active TSS

    if stitch_window == "":
        print("DETERMINING OPTIMUM STITCHING PARAMTER")
        opt_collection = copy.deepcopy(reference_collection)
        stitch_window = optimize_stitching(
            opt_collection, name, out_folder, step_size=500
        )
    print("USING A STITCHING PARAMETER OF {}".format(stitch_window))
    stitched_collection = reference_collection.stitch_collection(stitch_window, "both")

    if remove_tss:
        # now replace any stitched region that overlap 2 distinct genes
        # with the original loci that were there
        fixed_loci = []
        tss_loci = []
        for gene_id in list(start_dict.keys()):
            tss_loci.append(utils.make_tss_locus(gene_id, start_dict, 50, 50))

        # this turns the tss_loci list into a LocusCollection
        # 50 is the internal parameter for LocusCollection and doesn't really matter
        tss_collection = utils.LocusCollection(tss_loci, 50)
        remove_ticker = 0
        original_ticker = 0
        for stitched_locus in stitched_collection.get_loci():
            overlapping_tss_loci = tss_collection.get_overlap(stitched_locus, "both")
            tss_names = [
                start_dict[tss_locus.id]["name"] for tss_locus in overlapping_tss_loci
            ]
            tss_names = utils.uniquify(tss_names)
            if len(tss_names) > 2:

                # stitched_collection.remove(stitched_locus)
                original_loci = reference_collection.get_overlap(stitched_locus, "both")
                original_ticker += len(original_loci)
                fixed_loci += original_loci
                debug_output.append(
                    [stitched_locus.__str__(), stitched_locus.id, "MULTIPLE_TSS"]
                )
                remove_ticker += 1
            else:
                fixed_loci.append(stitched_locus)

        print(
            "REMOVED {} STITCHED LOCI BECAUSE THEY OVERLAPPED MULTIPLE TSSs".format(
                str(remove_ticker)
            )
        )
        print("ADDED BACK {} ORIGINAL LOCI".format(str(original_ticker)))
        fixed_collection = utils.LocusCollection(fixed_loci, 50)
        return fixed_collection, debug_output, stitch_window
    else:
        return stitched_collection, debug_output, stitch_window


# ==================================================================
# =====================REGION LINKING MAPPING=======================
# ==================================================================


def map_collection(
    stitched_collection,
    reference_collection,
    bam_file_list,
    mapped_folder,
    output,
    ref_name,
):
    """Makes a table of factor density in a stitched locus.

    Rank table by number of loci stitched together.

    """
    print("FORMATTING TABLE")
    loci = list(stitched_collection.get_loci())

    locus_table = [
        ["REGION_ID", "CHROM", "START", "STOP", "NUM_LOCI", "CONSTITUENT_SIZE"]
    ]

    loci_len_list = []

    # strip out any that are in chrY
    for locus in loci:
        if locus.chr == "chrY":
            loci.remove(locus)

    for locus in loci:
        # numLociList.append(int(stitchLocus.id.split('_')[1]))
        loci_len_list.append(locus.len())
        # numOrder = order(numLociList,decreasing=True)
    len_order = utils.order(loci_len_list, decreasing=True)
    ticker = 0
    for i in len_order:
        ticker += 1
        if ticker % 1000 == 0:
            print(ticker)
        locus = loci[i]

        # First get the size of the enriched regions within the stitched locus
        ref_enrich_size = 0
        ref_overlapping_loci = reference_collection.get_overlap(locus, "both")
        for ref_locus in ref_overlapping_loci:
            ref_enrich_size += ref_locus.len()

        try:
            stitch_count = int(locus.id.split("_")[0])
        except ValueError:
            stitch_count = 1
        coords = [int(x) for x in locus.coords()]

        locus_table.append(
            [
                locus.id,
                locus.chr,
                min(coords),
                max(coords),
                stitch_count,
                ref_enrich_size,
            ]
        )

    print("GETTING MAPPED DATA")
    print("USING A bam_file LIST:")
    print(bam_file_list)
    for bam_file in bam_file_list:

        bam_file_name = os.path.basename(bam_file)

        print("GETTING MAPPING DATA FOR  {}".format(bam_file))
        # assumes standard convention for naming enriched region gffs

        # opening up the mapped GFF
        mapped_gff_file = os.path.join(
            mapped_folder, "{}_{}_MAPPED".format(ref_name, bam_file_name), "matrix.txt"
        )
        print("OPENING {}".format(mapped_gff_file))

        mapped_gff = utils.parse_table(mapped_gff_file, "\t")

        signal_dict = defaultdict(float)
        print("MAKING SIGNAL DICT FOR {}".format(bam_file))
        mapped_loci = []
        for line in mapped_gff[1:]:

            chrom = line[1].split("(")[0]
            start = int(line[1].split(":")[-1].split("-")[0])
            end = int(line[1].split(":")[-1].split("-")[1])
            mapped_loci.append(utils.Locus(chrom, start, end, ".", line[0]))
            try:
                signal_dict[line[0]] = float(line[2]) * (abs(end - start))
            except ValueError:
                print("WARNING NO SIGNAL FOR LINE:")
                print(line)
                continue

        mapped_collection = utils.LocusCollection(mapped_loci, 500)
        locus_table[0].append(bam_file_name)

        for i in range(1, len(locus_table)):
            signal = 0.0
            line = locus_table[i]
            line_locus = utils.Locus(line[1], line[2], line[3], ".")
            overlapping_regions = mapped_collection.get_overlap(
                line_locus, sense="both"
            )
            for region in overlapping_regions:
                signal += signal_dict[region.id]
            locus_table[i].append(signal)

    utils.unparse_table(locus_table, output, "\t")
