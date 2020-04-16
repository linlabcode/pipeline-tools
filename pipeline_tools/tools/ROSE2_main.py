#!/usr/bin/env python

"""
PROGRAM TO STITCH TOGETHER REGIONS TO FORM ENHANCERS, MAP READ DENSITY TO STITCHED REGIONS,
AND RANK ENHANCERS BY READ DENSITY TO DISCOVER SUPER-ENHANCERS
"""

import argparse
import os
import sys
import time
from shutil import copyfile


from pipeline_tools.definitions import ROOT_DIR
from pipeline_tools.utils import rose2_utils, utils


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
        help="Enter a .gff or .bed file of binding sites used to make enhancers",
    )
    parser.add_argument(
        "-r",
        "--rankby",
        dest="rankby",
        required=True,
        help="bam_file to rank enhancer by",
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
        "-b",
        "--bams",
        dest="bams",
        required=False,
        help="Enter a comma separated list of additional bam files to map to",
    )
    parser.add_argument(
        "-c",
        "--control",
        dest="control",
        required=False,
        help="bam_file to rank enhancer by",
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
    mapped_folder = utils.format_folder(os.path.join(out_folder, "mapped_gff"), True)

    # GETTING INPUT FILE
    if args.input.split(".")[-1] == "bed":
        # CONVERTING A BED TO GFF
        input_gff_name = args.input.split("/")[-1][0:-4]
        input_gff_file = os.path.join(gff_folder, "{}.gff".format(input_gff_name))
        utils.bed_to_gff(args.input, input_gff_file)
    elif args.input.split(".")[-1] == "gff":
        # COPY THE INPUT GFF TO THE GFF FOLDER
        input_gff_file = args.input
        copyfile(
            input_gff_file, os.path.join(gff_folder, os.path.basename(input_gff_file))
        )

    else:
        print(
            "WARNING: INPUT FILE DOES NOT END IN .gff or .bed. ASSUMING .gff FILE FORMAT"
        )
        # COPY THE INPUT GFF TO THE GFF FOLDER
        input_gff_file = args.input
        copyfile(
            input_gff_file, os.path.join(gff_folder, os.path.basename(input_gff_file))
        )

    # GETTING THE LIST OF bam_fileS TO PROCESS
    if args.control:
        bam_file_list = [args.rankby, args.control]

    else:
        bam_file_list = [args.rankby]

    if args.bams:
        bam_file_list += args.bams.split(",")
        # bam_file_list = utils.uniquify(bam_file_list) # makes sad when you have the same control
        # bam over and over again
    # optional args

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

    # GETTING THE BOUND REGION FILE USED TO DEFINE ENHANCERS
    print("USING {} AS THE INPUT GFF".format(input_gff_file))
    input_name = os.path.basename(input_gff_file).split(".")[0]

    # GETTING THE GENOME
    genome = args.genome
    print("USING {} AS THE GENOME".format(genome))

    annot_file = rose2_utils.genome_dict[genome.upper()]

    # GET CHROMS FOUND IN THE BAMS
    print("GETTING CHROMS IN bam_fileS")
    bam_chrom_list = rose2_utils.get_bam_chrom_list(bam_file_list)
    print("USING THE FOLLOWING CHROMS")
    print(bam_chrom_list)

    # LOADING IN THE GFF AND FILTERING BY CHROM
    print("LOADING AND FILTERING THE GFF")
    input_gff = rose2_utils.filter_gff(input_gff_file, bam_chrom_list)
    # LOADING IN THE BOUND REGION REFERENCE COLLECTION
    print("LOADING IN GFF REGIONS")
    reference_collection = utils.gff_to_locus_collection(input_gff)
    print("STARTING WITH {} INPUT REGIONS".format(len(reference_collection)))
    print("CHECKING REFERENCE COLLECTION:")
    rose2_utils.check_ref_collection(reference_collection)

    # MASKING REFERENCE COLLECTION
    # see if there's a mask
    if args.mask:
        mask_file = args.mask
        print("USING MASK FILE {}".format(mask_file))
        # if it's a bed file
        if mask_file.split(".")[-1].upper() == "BED":
            mask_gff = utils.bed_to_gff(mask_file)
        elif mask_file.split(".")[-1].upper() == "GFF":
            mask_gff = utils.parse_table(mask_file, "\t")
        else:
            print("MASK MUST BE A .gff or .bed FILE")

        mask_collection = utils.gff_to_locus_collection(mask_gff)
        print("LOADING {} MASK REGIONS".format(str(len(mask_collection))))
        # now mask the reference loci
        reference_loci = reference_collection.get_loci()
        filtered_loci = [
            locus
            for locus in reference_loci
            if len(mask_collection.get_overlap(locus, "both")) == 0
        ]
        print(
            "FILTERED OUT {} LOCI THAT WERE MASKED IN {}".format(
                str(len(reference_loci) - len(filtered_loci)), mask_file
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
    # making sure start/stop ordering are correct
    for i in range(len(stitched_gff)):

        line = stitched_gff[i]
        start = int(line[3])
        stop = int(line[4])
        if start > stop:
            line[3] = stop
            line[4] = start

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

    print("CALLING AND PLOTTING SUPER-ENHANCERS")

    if args.control:
        control_name = os.path.basename(args.control)
    else:
        control_name = "NONE"
    cmd = "Rscript {} {} {} {} {}".format(
        os.path.join(ROOT_DIR, "scripts", "ROSE2_callSuper.R"),
        out_folder + "/",  # TODO: fix R script so it does not require '/'
        output_file1,
        input_name,
        control_name,
    )
    print(cmd)

    os.system(cmd)

    # calling the gene mapper
    time.sleep(20)
    super_table_file = "{}_SuperEnhancers.table.txt".format(input_name)
    if args.control:
        cmd = "ROSE2_geneMapper -g {} -r {} -c {} -i {}".format(
            genome,
            args.rankby,
            args.control,
            os.path.join(out_folder, super_table_file),
        )
    else:
        cmd = "ROSE2_geneMapper -g {} -r {} -i {}".format(
            genome, args.rankby, os.path.join(out_folder, super_table_file)
        )
    os.system(cmd)

    stretch_table_file = "{}_StretchEnhancers.table.txt".format(input_name)
    if args.control:
        cmd = "ROSE2_geneMapper -g {} -r {} -c {} -i {}".format(
            genome,
            args.rankby,
            args.control,
            os.path.join(out_folder, stretch_table_file),
        )
    else:
        cmd = "ROSE2_geneMapper -g {} -r {} -i {}".format(
            genome, args.rankby, os.path.join(out_folder, stretch_table_file)
        )
    os.system(cmd)

    superstretch_table_file = "{}_SuperStretchEnhancers.table.txt".format(input_name)
    if args.control:
        cmd = "ROSE2_geneMapper -g {} -r {} -c {} -i {}".format(
            genome,
            args.rankby,
            args.control,
            os.path.join(out_folder, superstretch_table_file),
        )
    else:
        cmd = "ROSE2_geneMapper -g {} -r {} -i {}".format(
            genome, args.rankby, os.path.join(out_folder, superstretch_table_file)
        )
    os.system(cmd)


if __name__ == "__main__":
    main()
