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

import os
import sys

import numpy
from collections import defaultdict
from pipeline_tools.utils import utils


# ============================================================================
# ============================MAPPING FUNCTIONS===============================
# ============================================================================

# dataFile -> bam_list
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
