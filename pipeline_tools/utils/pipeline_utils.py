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

import datetime
import os
import random
import sys

import numpy
from collections import defaultdict
from pipeline_tools.utils import utils


# ============================================================================
# ============================MAPPING FUNCTIONS===============================
# ============================================================================

# data_file -> bam_list
def map_regions(
    bam_list,
    gff_list,
    mapped_folder,
    signal_folder,
    names_list=[],
    median_norm=False,
    output="",
    extend_reads_to=200,
):
    """Making a normalized binding signal table at all regions."""
    # set up a list to return all signal tables made
    signal_table_list = []
    # since each bam has different read lengths, important to carefully normalize quantification
    if not names_list:
        names_list = [os.path.basename(bam).split(".")[0] for bam in bam_list]

    print(names_list)
    for name, bam_file in sorted(zip(names_list, bam_list)):
        bam = utils.Bam(bam_file)
        read_length = bam.get_read_lengths()[0]
        if int(extend_reads_to) < read_length:
            print(
                "Error: desired overall read extension {} is less than read length {}".format(
                    extend_reads_to, read_length,
                )
            )
            sys.exit()

        bam_extension = int(extend_reads_to) - read_length
        print("For dataset {} using an extension of {}".format(name, bam_extension))
        map_bams_batch(
            [bam_file],
            gff_list,
            mapped_folder,
            overwrite=False,
            names_list=[name],
            extension=bam_extension,
            rpm=True,
        )

    # want a signal table of all datasets to each gff
    print("Writing signal tables for each gff:")
    for gff_file in gff_list:
        gff_name = os.path.basename(gff_file).split(".")[0]
        if not output:
            signal_table_path = os.path.join(
                signal_folder, "{}_SIGNAL.txt".format(gff_name)
            )
        else:
            signal_table_path = output
        print(signal_table_path)
        make_signal_table(
            names_list, gff_file, mapped_folder, median_norm, output=signal_table_path,
        )
        signal_table_list.append(signal_table_path)

    return signal_table_list


def map_bams_batch(
    bam_list,
    gff_list,
    mapped_folder,
    overwrite=False,
    names_list=[],
    extension=200,
    rpm=True,
):
    """For each gff maps all of the data and writes to a specific folder named after the gff.

    Can map either by cell type or by a specific name list.
    Uses bamliquidator_batch.

    """
    if not names_list:
        names_list = [os.path.basename(bam).split(".")[0] for bam in bam_list]

    for gff_file in gff_list:
        # check to make sure gff exists
        try:
            open(gff_file, "r").close()
        except IOError:
            print("ERROR: GFF FILE {} DOES NOT EXIST".format(gff_file))
            sys.exit()

        gff_name = os.path.basename(gff_file).split(".")[0]

        # see if the parent directory exists, if not make it
        mapped_folder = utils.format_folder(mapped_folder, True)
        outdir_root = utils.format_folder(os.path.join(mapped_folder, gff_name), True)

        for name, bam in zip(names_list, bam_list):
            print("mapping {} to {}".format(name, gff_file))

            # filter based on celltype
            # output for the bamliquidator command
            outdir = utils.format_folder(os.path.join(outdir_root, name), True)
            out_matrix_file = os.path.join(outdir, "matrix.txt")

            if overwrite:
                map_cmd = "bamliquidator_batch --sense . -e {} --match_bamToGFF -r {} -o {} {}".format(
                    extension, gff_file, outdir, bam,
                )
                print(map_cmd)
                os.system(map_cmd)

            else:
                try:
                    print("checking for outfile {}".format(out_matrix_file))
                    open(out_matrix_file, "r").close()
                    print("File {} Already Exists, not mapping".format(out_matrix_file))
                except IOError:
                    map_cmd = "bamliquidator_batch --sense . -e {} --match_bamToGFF -r {} -o {} {}".format(
                        extension, gff_file, outdir, bam,
                    )
                    print(map_cmd)
                    os.system(map_cmd)

    # now initiate another giant loop to check for output and rename it
    for gff_file in gff_list:
        # check to make sure gff exists
        try:
            open(gff_file, "r").close()
        except IOError:
            print("ERROR: GFF FILE {} DOES NOT EXIST".format(gff_file))
            sys.exit()

        gff_name = os.path.basename(gff_file).split(".")[0]

        # see if the parent directory exists, if not make it
        mapped_folder = utils.format_folder(mapped_folder, True)
        # the first outdir of the mapping
        outdir_root = utils.format_folder(os.path.join(mapped_folder, gff_name), True)

        for name in names_list:
            print("Checking output of {} mapping to {}".format(name, gff_file))

            outdir = utils.format_folder(os.path.join(outdir_root, name), True)
            matrix_file = os.path.join(outdir, "matrix.txt")

            # what we want the eventual outfile to look like
            out_matrix_file = os.path.join(
                outdir_root, "{}_{}.txt".format(gff_name, name)
            )

            # now make sure the matrix file exists
            try:
                open(out_matrix_file, "r").close()
            except IOError:
                if utils.check_output(matrix_file, 0.1, 2):
                    print(
                        "Renaming output {} as {}".format(matrix_file, out_matrix_file)
                    )
                    os.rename(matrix_file, out_matrix_file)
                else:
                    print(
                        "ERROR: No output found for {} mapping to {}".format(
                            name, gff_file
                        )
                    )


def make_signal_table(
    names_list, gff_file, mapped_folder, median_norm=False, output=""
):
    """For each sample, make a dictionary keyed by locus ID."""
    signal_dict = {}
    for name in names_list:
        signal_dict[name] = defaultdict(float)

    # now start filling in the signal dict
    gff_name = os.path.basename(gff_file).split(".")[0]
    print(gff_name)
    for name in names_list:
        print("MAKING SIGNAL DICT FOR %s" % (name))

        # try opening the batch mapping output first
        mapped_file = os.path.join(
            mapped_folder, gff_name, "{}_{}.txt".format(gff_name, name)
        )
        if utils.check_output(mapped_file, 0.02, 0.02):
            print("FOUND MAPPED FILE FOR {} AT {}".format(name, mapped_file))
        else:
            mapped_file = os.path.join(
                mapped_folder, gff_name, "{}_{}.txt".format(gff_name, name),
            )

        if utils.check_output(mapped_file, 0.02, 0.02):
            print("FOUND MAPPED FILE FOR {} AT {}".format(name, mapped_file))
        else:
            print("ERROR NO MAPPED FILE FOUND FOR {}".format(name))
            sys.exit()

        mapped_table = utils.parse_table(mapped_file, "\t")
        if median_norm:
            median_signal = numpy.median([float(line[2]) for line in mapped_table[1:]])
        else:
            median_signal = 1

        for line in mapped_table[1:]:
            signal_dict[name][line[1]] = float(line[2]) / median_signal

    # now make the signal table
    signal_table = []
    header = ["GENE_ID", "locusLine"] + names_list
    signal_table.append(header)

    for line in mapped_table[1:]:
        locus_id = line[1]
        sig_line = line[0:2] + [signal_dict[name][locus_id] for name in names_list]
        signal_table.append(sig_line)

    if not output:
        return signal_table
    else:
        utils.unparse_table(signal_table, output, "\t")
        return signal_table


# =================================================================
# ===================DATA TABLE FUNCTIONS==========================
# =================================================================


def load_data_table(data_file):
    """Load the master data table."""

    if isinstance(data_file, str):
        data_table = utils.parse_table(data_file, "\t")
    else:
        data_table = list(data_file)
    # first check to make sure the table is formatted correctly
    for line in data_table:
        # print(line)
        if len(line) != 9:
            print("this line did not pass")
            print(line)
            data_table = format_data_table(data_file)
            break

    data_dict = defaultdict(dict)
    for line in data_table[1:]:
        data_dict[line[3]]["folder"] = utils.format_folder(line[0], False)
        data_dict[line[3]]["uniqueID"] = line[1]
        data_dict[line[3]]["genome"] = line[2].upper()
        genome = line[2]

        data_dict[line[3]]["sam"] = "".join([line[0], line[1], ".", genome, ".bwt.sam"])
        data_dict[line[3]]["ylf"] = "".join([line[0], line[1], ".", genome, ".bwt.ylf"])
        data_dict[line[3]]["enriched"] = line[5]
        data_dict[line[3]]["background"] = line[4]
        data_dict[line[3]]["enrichedMacs"] = line[6]
        color_string = line[7].replace('"', "")
        data_dict[line[3]]["color"] = color_string
        data_dict[line[3]]["fastq"] = line[8]

        # figure out which bam convention we are using
        # default will be new convention
        # look in the bam_folder for all bams that might fit the bill
        bam_folder = str(line[0])
        bam_file_list = [
            x for x in os.listdir(bam_folder) if len(x) > 0 and x[0] != "."
        ]

        bam_file_candidates = [
            x
            for x in bam_file_list
            if x.count(line[1]) == 1
            and x.split(".")[-1] == "bam"
            and x.count("bai") == 0
        ]
        if not bam_file_candidates:
            print(
                "UNABLE TO FIND A BAM FILE IN {} WITH UNIQUE ID {}".format(
                    bam_folder, line[1]
                )
            )
            full_bam_path = ""
        elif len(bam_file_candidates) > 1:
            print(
                "MUTLIPLE BAM FILES IN {} WITH UNIQUE ID {}. NO BAM ASISGNED".format(
                    bam_folder, line[1]
                )
            )
            print(bam_file_candidates)
            full_bam_path = ""
        else:
            bam_file = bam_file_candidates[0]
            full_bam_path = os.path.abspath(os.path.join(bam_folder, bam_file))
            full_bai_path = full_bam_path + ".bai"

        if full_bam_path:
            try:
                open(full_bam_path, "r").close()
            except (IOError, FileNotFoundError):
                print("ERROR: BAM FILE {} DOES NOT EXIST".format(full_bam_path))
                full_bam_path = ""
            try:
                open(full_bai_path, "r").close()
            except (IOError, FileNotFoundError):
                print(
                    "ERROR: BAM FILE {} DOES NOT HAVE BAI INDEX".format(full_bam_path)
                )
                full_bam_path = ""

        data_dict[line[3]]["bam"] = full_bam_path

    return data_dict


def format_data_table(data_file):
    """Formats the data_file and rewrite.

    First 3 columns are required for every line. If they aren't there the line is deleted.

    """
    print("reformatting data table")

    data_table = utils.parse_table(data_file, "\t")

    new_data_table = [
        [
            "FILE_PATH",
            "UNIQUE_ID",
            "GENOME",
            "NAME",
            "BACKGROUND",
            "ENRICHED_REGION",
            "ENRICHED_MACS",
            "COLOR",
            "FASTQ_FILE",
        ]
    ]
    # first check to make sure the table is formatted correctly
    for line in data_table[1:]:
        if len(line) < 3:
            continue
        # this spots header lines that may be out of place
        if line[0] == "FILE_PATH":
            continue
        # check if it at least has the first 3 columns filled in
        if len(line[0]) == 0 or len(line[1]) == 0 or len(line[2]) == 0:
            print("ERROR required fields missing in line")
            print(line)
        # if the first three are filled in, check to make sure there are 8 columns
        else:
            if len(line) > 3 and len(line) < 9:
                new_line = line + (8 - len(line)) * [""] + ["NA"]
                new_data_table.append(new_line)
            elif len(line) >= 9:
                new_line = line[0:9]
                new_data_table.append(new_line)

    # lower case all of the genomes
    # make the color 0,0,0 for blank lines and strip out any " marks
    for i in range(1, len(new_data_table)):
        new_data_table[i][2] = new_data_table[i][2].lower()
        color = new_data_table[i][7]
        if len(color) == 0:
            new_data_table[i][7] = "0,0,0"
    utils.unparse_table(new_data_table, data_file, "\t")

    return new_data_table


# ==========================================================================
# ==============================CALLING ROSE================================
# ==========================================================================


def call_rose2(
    data_file,
    macs_enriched_folder,
    parent_folder,
    names_list=[],
    extra_map=[],
    input_file="",
    tss=2500,
    stitch="",
    bash_file_name="",
    mask="",
    use_background=True,
):
    """Call rose w/ standard parameters."""
    data_dict = load_data_table(data_file)
    if len(names_list) == 0:
        names_list = list(data_dict.keys())

    # a time_stamp to name this pipeline batch of files
    time_stamp = datetime.datetime.now().strftime("%Y%m%d%H%M%S")
    # a random integer ticker to help name files
    rand_ticker = random.randint(0, 10000)

    utils.format_folder(parent_folder, True)

    if len(bash_file_name) == 0:
        bash_file_name = os.path.join(
            parent_folder, "rose_{}_{}.sh".format(time_stamp, rand_ticker)
        )
    with open(bash_file_name, "w") as bash_file:
        map_string = [data_dict[name]["bam"] for name in extra_map]
        map_string = ",".join(map_string)

        for name in names_list:
            genome = data_dict[name]["genome"]
            bam_file = data_dict[name]["bam"]

            background_name = data_dict[name]["background"]
            if use_background and background_name in data_dict:
                background_bam_file = data_dict[background_name]["bam"]
                has_background = True
            else:
                has_background = False

            if len(input_file) == 0:
                macs_file = os.path.join(
                    macs_enriched_folder, data_dict[name]["enrichedMacs"]
                )
            else:
                macs_file = input_file
            output_folder = os.path.join(parent_folder, "{}_ROSE".format(name))
            print(name)
            bash_file.write("#running ROSE2 on {}\n".format(name))
            rose_cmd = "ROSE2 -g {} -i {} -r {} -o {} -t {}".format(
                genome, macs_file, bam_file, output_folder, tss
            )

            if len(str(stitch)) > 0:
                rose_cmd += " -s {}".format(stitch)
            if has_background:
                rose_cmd += " -c {}".format(background_bam_file)
            if len(map_string) > 0:
                rose_cmd += " -b {}".format(map_string)
            if len(mask) > 0:
                rose_cmd += " --mask {}".format(mask)

            bash_file.write(rose_cmd)
            bash_file.write("\n\n")

    print("Wrote rose commands to {}".format(bash_file_name))

    return bash_file_name


# ==============================================================================
# ==============================ENHANCER FUNCTIONS==============================
# ==============================================================================


def get_file(file_string, file_list, parent_folder):
    """Returns full path of file from file_list containing the file_string.

    Return an error if multiple files match.

    """
    if not utils.format_folder(parent_folder, False):
        print("ERROR: Folder {} does not exist".format(parent_folder))
        sys.exit()
    parent_folder = utils.format_folder(parent_folder, False)
    match_files = [
        file_name for file_name in file_list if file_name.count(file_string) == 1
    ]
    if not match_files:
        print(
            "ERROR: No files found in {} with {} in title".format(
                parent_folder, file_string
            )
        )
        sys.exit()
    if len(match_files) > 1:
        print(
            "ERROR: Multiple files found in {} with {} in title".format(
                parent_folder, file_string
            )
        )
        sys.exit()
    match_file_path = os.path.join(parent_folder, match_files[0])

    return match_file_path
