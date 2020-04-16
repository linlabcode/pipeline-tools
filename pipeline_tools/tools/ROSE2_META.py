#!/usr/bin/env python

"""
PROGRAM TO STITCH TOGETHER REGIONS TO FORM ENHANCERS, MAP READ DENSITY TO STITCHED REGIONS,
AND RANK ENHANCERS BY READ DENSITY TO DISCOVER SUPER-ENHANCERS
"""

import argparse
import os
import sys
from shutil import copyfile

import numpy

from pipeline_tools.definitions import ROOT_DIR
from pipeline_tools.utils import rose2_utils, utils

# ==================================================================
# ====================COLLAPSING REGION MAP=========================
# ==================================================================


def collapse_region_map(region_map_file, name="", control_bams=False):
    """Take a region_map file and collapse signal into a single column.

    Also fix any stupid start/stop sorting issues. Need to take into account whether or not
    controls were used.

    """
    region_map = utils.parse_table(region_map_file, "\t")

    for n, line in enumerate(region_map):
        if n == 0:
            # new header
            if len(name) == 0:
                name = "MERGED_SIGNAL"
            region_map[n] = line[0:6] + [name]

        else:
            new_line = list(line[0:6])
            if control_bams:
                signal_line = [float(x) for x in line[6:]]
                rankby_indexes = range(0, len(signal_line) // 2, 1)
                control_indexes = range(len(signal_line) // 2, len(signal_line), 1)
                meta_vector = []
                for i, j in zip(rankby_indexes, control_indexes):
                    # min signal is 0
                    meta_vector.append(max(0, signal_line[i] - signal_line[j]))
                meta_signal = numpy.mean(meta_vector)
            else:
                meta_signal = numpy.mean([float(x) for x in line[6:]])
            region_map[n] = new_line + [meta_signal]

    output_file = region_map_file.replace("REGION", "META")
    utils.unparse_table(region_map, output_file, "\t")
    return output_file


# ==================================================================
# =========================MAIN METHOD==============================
# ==================================================================
def main():
    """Main run call."""
    debug = False

    parser = argparse.ArgumentParser()
    # required flags
    parser.add_argument(
        "-i",
        "--i",
        dest="input",
        required=True,
        help=(
            "Enter a comma separated list of .gff or .bed file of binding sites used to make "
            "enhancers"
        ),
    )
    parser.add_argument(
        "-r",
        "--rankby",
        dest="rankby",
        required=True,
        help="Enter a comma separated list of bams to rank by",
    )
    parser.add_argument(
        "-o", "--out", dest="out", required=True, help="Enter an output folder"
    )
    parser.add_argument(
        "-g",
        "--genome",
        dest="genome",
        required=True,
        help="Enter the genome build (MM9,MM8,HG18,HG19)",
    )

    # optional flags
    parser.add_argument(
        "-n",
        "--name",
        dest="name",
        required=False,
        help="Provide a name for the analysis otherwise ROSE will guess",
    )
    parser.add_argument(
        "-c",
        "--control",
        dest="control",
        required=False,
        help=(
            "Enter a comma separated list of control bams. Can either provide a single control "
            "bam for all rankby bams, or provide a control bam for each individual bam"
        ),
    )
    parser.add_argument(
        "-s",
        "--stitch",
        dest="stitch",
        default="",
        help=(
            "Enter a max linking distance for stitching. Default will determine optimal stitching"
            " parameter"
        ),
    )
    parser.add_argument(
        "-t",
        "--tss",
        dest="tss",
        default=0,
        help="Enter a distance from TSS to exclude. 0 = no TSS exclusion",
    )

    parser.add_argument(
        "--mask",
        dest="mask",
        required=False,
        help="Mask a set of regions from analysis.  Provide a .bed or .gff of masking regions",
    )

    # RETRIEVING FLAGS
    args = parser.parse_args()

    # making the out folder if it doesn't exist
    out_folder = utils.format_folder(args.out, True)

    # figuring out folder schema
    gff_folder = utils.format_folder(os.path.join(out_folder, "gff"), True)
    mapped_folder = utils.format_folder(os.path.join(out_folder, "mappedGFF"), True)

    # GETTING INPUT FILE(s)
    input_list = [
        input_file for input_file in args.input.split(",") if len(input_file) > 1
    ]

    # converting all input files into GFFs and moving into the GFF folder
    input_gf_list = []
    for input_file in input_list:
        # GETTING INPUT FILE
        if args.input.split(".")[-1] == "bed":
            # CONVERTING A BED TO GFF
            input_gff_name = os.path.basename(args.input)[0:-4]
            input_gff_file = os.path.join(gff_folder, "{}.gff".format(input_gff_name))
            utils.bed_to_gff(args.input, input_gff_file)
        elif args.input.split(".")[-1] == "gff":
            # COPY THE INPUT GFF TO THE GFF FOLDER
            input_gff_file = args.input
            copyfile(
                input_gff_file,
                os.path.join(gff_folder, os.path.basename(input_gff_file)),
            )
        else:
            print(
                "WARNING: INPUT FILE DOES NOT END IN .gff or .bed. ASSUMING .gff FILE FORMAT"
            )
            # COPY THE INPUT GFF TO THE GFF FOLDER
            input_gff_file = args.input
            copyfile(
                input_gff_file,
                os.path.join(gff_folder, os.path.basename(input_gff_file)),
            )

        input_gf_list.append(input_gff_file)

    # GETTING THE LIST OF bam_fileS TO PROCESS
    # either same number of bams for rankby and control
    # or only 1 control #or none!
    # bamlist should be all rankby bams followed by control bams

    bam_file_list = []
    if args.control:
        control_bam_list = [bam for bam in args.control.split(",") if len(bam) > 0]
        rankby_bam_list = [bam for bam in args.rankby.split(",") if len(bam) > 0]

        if len(control_bam_list) == len(rankby_bam_list):
            # case where an equal number of backgrounds are given
            bam_file_list = rankby_bam_list + control_bam_list
        elif len(control_bam_list) == 1:
            # case where a universal background is applied
            bam_file_list = rankby_bam_list + control_bam_list * len(rankby_bam_list)
        else:
            print(
                "ERROR: EITHER PROVIDE A SINGLE CONTROL BAM FOR ALL SAMPLES, OR ONE CONTROL BAM"
                " FOR EACH SAMPLE"
            )
            sys.exit()
    else:
        bam_file_list = [bam for bam in args.rankby.split(",") if len(bam) > 0]

    # Stitch parameter
    if args.stitch == "":
        stitch_window = ""
    else:
        stitch_window = int(args.stitch)

    # tss args
    tss_window = int(args.tss)
    if tss_window != 0:
        remove_tss = True
    else:
        remove_tss = False

    # GETTING THE GENOME
    genome = args.genome.upper()
    print("USING {} AS THE GENOME".format(genome))

    # GETTING THE CORRECT ANNOT FILE
    try:
        annot_file = rose2_utils.genome_dict[genome]
    except KeyError:
        print("ERROR: UNSUPPORTED GENOMES TYPE {}".format(genome))
        sys.exit()

    # FINDING THE ANALYSIS NAME
    if args.name:
        input_name = args.name
    else:
        input_name = os.path.basename(input_gf_list[0]).split(".")[0]
    print("USING {} AS THE ANALYSIS NAME".format(input_name))

    print("FORMATTING INPUT REGIONS")
    # MAKING THE RAW INPUT FILE FROM THE INPUT GFFs
    # use a simpler unique region naming system
    if len(input_gf_list) == 1:
        input_gff = utils.parse_table(input_gf_list[0], "\t")
    else:
        input_loci = []
        for gff_file in input_gf_list:
            print("\tprocessing {}".format(gff_file))
            gff = utils.parse_table(gff_file, "\t")
            gff_collection = utils.gff_to_locus_collection(gff, 50)
            input_loci += gff_collection.get_loci()

        input_collection = utils.LocusCollection(input_loci, 50)
        input_collection = (
            input_collection.stitch_collection()
        )  # stitches to produce unique regions

        input_gff = utils.locus_collection_to_gff(input_collection)

    formatted_gff = []
    # now number things appropriately
    for i, line in enumerate(input_gff):

        # use the coordinates to make a new id input_name_chr_sense_start_stop
        chrom = line[0]
        coords = [int(line[3]), int(line[4])]
        sense = line[6]

        line_id = "{}_{}".format(input_name, str(i + 1))  # 1 indexing

        new_line = [
            chrom,
            line_id,
            line_id,
            min(coords),
            max(coords),
            "",
            sense,
            "",
            line_id,
        ]
        formatted_gff.append(new_line)

    # name of the master input gff file
    master_gff_file = os.path.join(
        gff_folder, "{}_{}_ALL_-0_+0.gff".format(genome, input_name)
    )
    utils.unparse_table(formatted_gff, master_gff_file, "\t")

    print("USING {} AS THE INPUT GFF".format(master_gff_file))

    # GET CHROMS FOUND IN THE BAMS
    print("GETTING CHROMS IN bam_fileS")
    bam_chrom_list = rose2_utils.get_bam_chrom_list(bam_file_list)
    print("USING THE FOLLOWING CHROMS")
    print(bam_chrom_list)

    # LOADING IN THE GFF AND FILTERING BY CHROM
    print("LOADING AND FILTERING THE GFF")
    input_gff = rose2_utils.filter_gff(master_gff_file, bam_chrom_list)
    # LOADING IN THE BOUND REGION REFERENCE COLLECTION
    print("LOADING IN GFF REGIONS")
    reference_collection = utils.gff_to_locus_collection(input_gff)

    print("CHECKING REFERENCE COLLECTION:")
    rose2_utils.check_ref_collection(reference_collection)

    # MASKING REFERENCE COLLECTION
    # see if there's a mask
    if args.mask:
        mask_file = args.mask
        # if it's a bed file
        if mask_file.split(".")[-1].upper() == "BED":
            mask_gff = utils.bedToGFF(mask_file)
        elif mask_file.split(".")[-1].upper() == "GFF":
            mask_gff = utils.parse_table(mask_file, "\t")
        else:
            print("MASK MUST BE A .gff or .bed FILE")
            sys.exit()
        mask_collection = utils.gff_to_locus_collection(mask_gff)

        # now mask the reference loci
        reference_loci = reference_collection.get_loci()
        filtered_loci = [
            locus
            for locus in reference_loci
            if len(mask_collection.get_overlap(locus, "both")) == 0
        ]
        print(
            "FILTERED OUT {} LOCI THAT WERE MASKED IN {}".format(
                len(reference_loci) - len(filtered_loci), mask_file
            )
        )
        reference_collection = utils.LocusCollection(filtered_loci, 50)

    # NOW STITCH REGIONS
    print("STITCHING REGIONS TOGETHER")
    stitched_collection, debug_output, stitch_window = rose2_utils.region_stitching(
        reference_collection,
        input_name,
        out_folder,
        stitch_window,
        tss_window,
        annot_file,
        remove_tss,
    )

    # NOW MAKE A STITCHED COLLECTION GFF
    print("MAKING GFF FROM STITCHED COLLECTION")
    stitched_gff = utils.locus_collection_to_gff(stitched_collection)

    print(stitch_window)
    print(type(stitch_window))
    if not remove_tss:
        stitched_gff_file = os.path.join(
            gff_folder,
            "{}_{}KB_STITCHED.gff".format(input_name, str(stitch_window // 1000)),
        )
        stitched_gff_name = "{}_{}KB_STITCHED".format(
            input_name, str(stitch_window // 1000)
        )
        debug_out_file = os.path.join(
            gff_folder,
            "{}_{}KB_STITCHED.debug".format(input_name, str(stitch_window // 1000)),
        )
    else:
        stitched_gff_file = os.path.join(
            gff_folder,
            "{}_{}KB_STITCHED_TSS_DISTAL.gff".format(
                input_name, str(stitch_window // 1000)
            ),
        )
        stitched_gff_name = "{}_{}KB_STITCHED_TSS_DISTAL".format(
            input_name, str(stitch_window // 1000)
        )
        debug_out_file = os.path.join(
            gff_folder,
            "{}_{}KB_STITCHED_TSS_DISTAL.debug".format(
                input_name, str(stitch_window // 1000)
            ),
        )

    # WRITING DEBUG OUTPUT TO DISK

    if debug:
        print("WRITING DEBUG OUTPUT TO DISK AS {}".format(debug_out_file))
        utils.unparse_table(debug_output, debug_out_file, "\t")

    # WRITE THE GFF TO DISK
    print("WRITING STITCHED GFF TO DISK AS {}".format(stitched_gff_file))
    utils.unparse_table(stitched_gff, stitched_gff_file, "\t")

    # SETTING UP THE OVERALL OUTPUT FILE
    output_file1 = os.path.join(
        out_folder, "{}_ENHANCER_REGION_MAP.txt".format(stitched_gff_name)
    )
    print("OUTPUT WILL BE WRITTEN TO  {}".format(output_file1))

    # MAPPING TO THE NON STITCHED (ORIGINAL GFF)
    # MAPPING TO THE STITCHED GFF

    bam_file_list_unique = list(bam_file_list)
    bam_file_list_unique = utils.uniquify(bam_file_list_unique)
    # prevent redundant mapping
    print("MAPPING TO THE FOLLOWING BAMS:")
    print(bam_file_list_unique)
    for bam_file in bam_file_list_unique:

        bam_file_name = os.path.basename(bam_file)

        # MAPPING TO THE STITCHED GFF
        mapped_out1_folder = os.path.join(
            mapped_folder, "{}_{}_MAPPED".format(stitched_gff_name, bam_file_name)
        )
        mapped_out1_file = os.path.join(
            mapped_folder,
            "{}_{}_MAPPED".format(stitched_gff_name, bam_file_name),
            "matrix.txt",
        )
        if utils.check_output(mapped_out1_file, 0.2, 0.2):
            print(
                "FOUND {} MAPPING DATA FOR BAM: {}".format(
                    stitched_gff_file, mapped_out1_file
                )
            )
        else:
            cmd1 = "bamliquidator_batch --sense . -e 200 --match_bamToGFF -r {} -o {} {}".format(
                stitched_gff_file, mapped_out1_folder, bam_file,
            )
            print(cmd1)

            os.system(cmd1)
            if utils.check_output(mapped_out1_file, 0.2, 5):
                print(
                    "SUCCESSFULLY MAPPED TO {} FROM BAM: {}".format(
                        stitched_gff_file, bam_file_name
                    )
                )
            else:
                print(
                    "ERROR: FAILED TO MAP {} FROM BAM: {}".format(
                        stitched_gff_file, bam_file_name
                    )
                )
                sys.exit()

    print("BAM MAPPING COMPLETED NOW MAPPING DATA TO REGIONS")
    # CALCULATE DENSITY BY REGION
    # NEED TO FIX THIS FUNCTION TO ACCOUNT FOR DIFFERENT OUTPUTS OF LIQUIDATOR
    rose2_utils.map_collection(
        stitched_collection,
        reference_collection,
        bam_file_list,
        mapped_folder,
        output_file1,
        ref_name=stitched_gff_name,
    )

    print("FINDING AVERAGE SIGNAL AMONGST BAMS")
    meta_output_file = collapse_region_map(
        output_file1, input_name + "_MERGED_SIGNAL", control_bams=args.control
    )

    # now try the merging

    print("CALLING AND PLOTTING SUPER-ENHANCERS")

    control_name = "NONE"
    cmd = "Rscript {} {} {} {} {}".format(
        os.path.join(ROOT_DIR, "scripts", "ROSE2_callSuper.R"),
        out_folder + "/",  # TODO: fix R script so it does not require '/'
        meta_output_file,
        input_name,
        control_name,
    )
    print(cmd)

    os.system(cmd)

    # calling the gene mapper
    print("CALLING GENE MAPPING")

    super_table_file = "{}_SuperEnhancers.table.txt".format(input_name)

    # for now don't use ranking bam to call top genes
    cmd = "ROSE2_geneMapper -g {} -i {} -f".format(
        genome, os.path.join(out_folder, super_table_file)
    )
    print(cmd)
    os.system(cmd)

    stretch_table_file = "{}_StretchEnhancers.table.txt".format(input_name)

    cmd = "ROSE2_geneMapper -g {} -i {} -f".format(
        genome, os.path.join(out_folder, stretch_table_file)
    )
    print(cmd)
    os.system(cmd)

    superstretch_table_file = "{}_SuperStretchEnhancers.table.txt".format(input_name)

    cmd = "ROSE2_geneMapper.py -g {} -i {} -f".format(
        genome, out_folder, superstretch_table_file
    )
    os.system(cmd)


if __name__ == "__main__":
    main()
