#!/bin/env python

import sys
import os
import os.path
import operator
import getopt
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "../"))
import naivetrack
import table


def variant_to_annovar(variant):
    """Convert a variant to the annovar format."""
    start = variant.start + 1
    end = variant.end + 1
    reference = variant.reference
    alternate = variant.alternate

    if len(alternate) > 1:
        ndel = alternate.count('-')
        if ndel > 0:
            reference = '0'
            alternate = '-'
            start += 1
            end += ndel
        else:
            reference = '-'
            alternate = alternate[1:]
            start += 1
            end += 1
    return table.Variant(variant.chrid, start, end, reference, alternate)


def naivetrack_to_annoinput(reader=None, writer=None):
    """Convert the naive tracks to the annovar input format."""

    # iterate over the data entries
    for entry in reader:

        # skip positions without data
        if len(entry.alleles) == 0:
            continue

        # the column data
        columns = {}

        # record amplicon information if present
        amplicons = {}
        if len(entry.alleles[0].tags) == 3:
            amplicondata = naivetrack.utils.aggregate_allele_information(
                entry, attributes=[
                    naivetrack.models.AlleleFactors["Allele"],
                    naivetrack.models.AlleleFactors["Sample"],
                    3])

            amplicons = {}
            for am in amplicondata.keys():
                amplicons[am] = {}
                for s in amplicondata[am].keys():
                    amplicons[am][s] = "%s" % (
                        ','.join(amplicondata[am][s].keys())
                    )

        # aggregate the data per allele, strand and sample
        per_allele = naivetrack.utils.aggregate_allele_information(
            entry, attributes=[
                naivetrack.models.AlleleFactors["Allele"],
                naivetrack.models.AlleleFactors["Sample"]])

        per_sample = naivetrack.utils.aggregate_allele_information(
            entry, attributes=[naivetrack.models.AlleleFactors["Sample"]])

        per_strand = naivetrack.utils.aggregate_allele_information(
            entry, attributes=[
                naivetrack.models.AlleleFactors["Allele"],
                naivetrack.models.AlleleFactors["Sample"],
                naivetrack.models.AlleleFactors["Strand"]])

        freq = naivetrack.utils.aggregate_divide(per_allele, per_sample)
        strandfreq = naivetrack.utils.aggregate_divide(per_strand, per_allele)

        for a in sorted(per_allele.keys()):
            variant = table.Variant(entry.position.chromosome,
                                    entry.position.coordinate,
                                    entry.position.coordinate,
                                    entry.position.sequence.upper(),
                                    a.upper())

            # prepare the data for writing
            columns[a] = {}
            samples = sorted(per_sample.keys())
            for s in samples:
                qualfreqkey = "01_Quality_frequency_[%s]" % s
                qualtotkey = "02_Quality_total_[%s]" % s
                depthfreqkey = "03_Depth_frequency_[%s]" % s
                depthaltkey = "04_Depth_alt_[%s]" % s
                depthrefkey = "05_Depth_ref_[%s]" % s
                depthtotkey = "06_Depth_total_[%s]" % s
                qualdepthkey = "07_Quality_Depth_[%s]" % s
                strandedness = "08_Strandedness_[%s]" % s
                ampkey = "09_Amplicons_[%s]" % s

                # no sample entries -> continue
                if s not in per_sample.keys():
                    continue

                # default values
                columns[a][qualtotkey] = per_sample[s]["__qual__"]
                columns[a][depthtotkey] = per_sample[s]["__depth__"]

                # no alternate allele -> 0.0
                try:
                    columns[a][qualfreqkey] = round(freq[a][s]["__qual__"], 4)
                    columns[a][depthfreqkey] = round(freq[a][s]["__depth__"], 4)
                except KeyError:
                    columns[a][qualfreqkey] = 0.0
                    columns[a][depthfreqkey] = 0.0
                except TypeError:
                    columns[a][qualfreqkey] = 0.0
                    columns[a][depthfreqkey] = 0.0

                # alternate allele found
                try:
                    columns[a][depthaltkey] = per_allele[a][s]["__depth__"]
                    columns[a][qualdepthkey] = round(operator.truediv(
                        per_allele[a][s]["__qual__"],
                        per_allele[a][s]["__depth__"]
                    ), 4)
                except KeyError:
                    columns[a][depthaltkey] = 0
                    columns[a][qualdepthkey] = 0.0
                except TypeError:
                    columns[a][qualdepthkey] = 0.0

                # reference allele found
                try:
                    columns[a][depthrefkey] = per_allele[variant.reference][s]["__depth__"]
                except KeyError:
                    columns[a][depthrefkey] = 0

                if amplicons is not None:
                    try:
                        columns[a][ampkey] = amplicons[a][s]
                    except KeyError:
                        columns[a][ampkey] = ""

                try:
                    columns[a][strandedness] = round(
                        strandfreq[a][s]['forward']["__qual__"], 4)
                except KeyError:
                    columns[a][strandedness] = 0.0
                except TypeError:
                    columns[a][strandedness] = 0.0

            if variant.alternate != variant.reference and len(variant.alternate) > 0:
                annovar = variant_to_annovar(variant)
                writer(annovar, columns[a])


if __name__ == "__main__":
    def usage(message="", error=None):
        """Print the usage information and messages."""
        sys.stderr.write(
"""
%(message)s

Called as: %(commandline)s

Usage:
    %(tool)s --input <inputfile.var> --output <outputfile.var>

Options:
    -i/--input                  [file]  The input var file, default stdin
    -o/--output                 [file]  The output annovar table, default stdout
    --minimum-quality           [int]   The minimum quality of a variant to include
    --maximum-quality           [int]   The maximum quality of a variant to include
    --minimum-quality-frequency [float] The minimum frequency in terms of quality of a variant to include
    --maximum-quality-frequency [float] The maximum frequency in terms of quality of a variant to include
    --minimum-depth             [int]   The minimum read depth of a variant to include
    --maximum-depth             [int]   The maximum read depth of a variant to include
    --minimum-depth-frequency   [float] The minimum frequency in terms of read depth of a variant to include
    --maximum-depth-frequency   [float] The maximum frequency in terms of read depth of a variant to include
    --remove-N                  []      Remove alternate alleles with N characters
    --no-reference              []      Remove all reference calls

""" % {
    "message": message,
    "commandline": " ".join(sys.argv),
    "tool": sys.argv[0]
    })
        if error is not None:
            sys.exit(int(error))

    def main():
        """Perform the main program loop."""
        # prepare data parsing
        chromosome_list = []
        provider = naivetrack.TrackParser(chromosome_list=chromosome_list)
        provider.set_maximum_chromosome_size(3000000000)

        # prepare the filter
        fltr = naivetrack.SampleFilter()
        fout = sys.stdout

        # parse the commandline options
        shortopt = "i:o:h"
        longopt = [
            "input=", "output=", "help",
            "minimum-quality=", "maximum-quality=",
            "minimum-quality-frequency=",
            "maximum-quality-frequency=",
            "minimum-depth=", "maximum-depth=",
            "remove-N", "no-reference",
            "minimum-depth-frequency=", "maximum-depth-frequency="
            ]
        try:
            opts, _ = getopt.getopt(sys.argv[1:], shortopt, longopt)
            for opt, ans in opts:
                if opt in ("-i", "--input"):
                    if ans != "-":
                        provider.register_handle(open(os.path.abspath(ans), 'r'))
                    else:
                        provider.register_handle(sys.stdin)
                elif opt in ("-o", "--output"):
                    fout = open(ans, 'w')
                elif opt == "--minimum-quality":
                    fltr.add_value_filter('__qual__', lowerbound=int(ans))
                elif opt == "--maximum-quality":
                    fltr.add_value_filter('__qual__', upperbound=int(ans))
                elif opt == "--minimum-quality-frequency":
                    fltr.add_frequency_filter('__qual__', lowerbound=float(ans))
                elif opt == "--maximum-quality-frequency":
                    fltr.add_frequency_filter('__qual__', upperbound=float(ans))
                elif opt == "--minimum-depth":
                    fltr.add_value_filter('__depth__', lowerbound=int(ans))
                elif opt == "--maximum-depth":
                    fltr.add_value_filter('__depth__', upperbound=int(ans))
                elif opt == "--minimum-depth-frequency":
                    fltr.add_frequency_filter('__depth__', lowerbound=float(ans))
                elif opt == "--maximum-depth-frequency":
                    fltr.add_frequency_filter('__depth__', upperbound=float(ans))
                elif opt == "--remove-N":
                    fltr.remove_N_alleles = True
                elif opt == "--no-reference":
                    fltr.remove_all_reference = True
                elif opt in ('-h', "--help"):
                    usage("Help was asked", error=0)
        except getopt.GetoptError as err:
            usage(str(err), error=101)

        # register filters
        if fltr.n_filters_defined() > 0:
            provider = naivetrack.apply_filters(provider=provider, filters=[fltr])

        # prepare the output writer
        outwriter = table.VariantWriter(fout, chromosome_list)

        # process the variants
        naivetrack_to_annoinput(provider, writer=outwriter)
        del outwriter

    # run the main program loop
    main()
