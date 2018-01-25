#!/bin/env python

import os
import os.path
import sys
import getopt
import itertools
import pysam


def count_amplicons(samin=None, quality=0, size=1000000):
    """Count amplicons in the input samfile."""
    assert samin is not None

    amplicons = {}
    amplicon_buffer = []
    for alignment in samin.fetch(until_eof=True):

        # don't consider unmapped alignments
        if alignment.is_unmapped:
            continue

        # don't consider alignments with a quality below the threshold
        if alignment.mapping_quality < quality:
            continue

        # get the amplicon
        tags = alignment.get_tags()
        amplicon = None
        for tag in tags:
            if tag[0] == "am":
                amplicon = tag[1]

        # don't consider reads without an amplicon
        if amplicon is None:
            continue

        # add the amplicon to the buffer
        amplicon_buffer.append(amplicon)

        # if the amplicon buffer is sufficiently large
        # add the amplicons to the return dict
        if len(amplicon_buffer) == size:
            amplicon_buffer.sort()
            for amplicon, itergrp in itertools.groupby(amplicon_buffer):
                count = 0
                for _ in itergrp:
                    count += 1
                if amplicon not in amplicons.keys():
                    amplicons[amplicon] = count
                else:
                    amplicons[amplicon] += count
            amplicon_buffer = []
            sys.stderr.write("%d Amplicons detected\n" % len(amplicons.keys()))

    # add the final amplicons to the return dict
    amplicon_buffer.sort()
    for amplicon, itergrp in itertools.groupby(amplicon_buffer):
        count = 0
        for _ in itergrp:
            count += 1
        if amplicon not in amplicons.keys():
            amplicons[amplicon] = count
        else:
            amplicons[amplicon] += count
    amplicon_buffer = []
    return amplicons


def design_from_bed(instream=sys.stdin):
    """Get the amplicon design from a BED file."""
    assert instream is not None

    design = []
    for line in instream:

        # remove the endline from the line
        line = line.rstrip()

        # filter bad lines
        if len(line) == 0:
            continue
        if line[0] == "#" or line.startswith("track"):
            continue

        # filter lines with fewer than 5 entries
        fields = line.split("\t")
        if len(fields) < 5:
            continue

        # format the amplicon
        amplicon = "%s:%s-%s(%s)" % (
            fields[0],
            fields[1],
            fields[2],
            fields[5]
        )
        design.append(amplicon)

    # return the amplicons in the design
    return design


def nimbus_count(samin=None, bedfile=sys.stdin, outstream=sys.stdout, quality=0):
    """Count the number of alignments per amplicon."""
    amplicons = count_amplicons(samin, quality)
    design = design_from_bed(bedfile)

    # print the amplicons in the design
    for amplicon in design:
        if amplicon in amplicons.keys():
            outstream.write("%s\t%d\n" % (amplicon, amplicons[amplicon]))
            del amplicons[amplicon]
        else:
            outstream.write("%s\t0\n" % (amplicon))

    # print the unknown amplicons
    unknown = 0
    for key, value in amplicons.items():
        unknown += value
        sys.stderr.write("UNKNOWN:%s\t%d\n" % (key, value))
    sys.stderr.write("%d Alignments for unknown amplicons detected\n" % unknown)


if __name__ == "__main__":

    def usage(message="", error=None):
        """Print the usage information."""
        sys.stdout.write(
"""
%(message)s

Called as: %(commandline)s

Usage:
    %(tool)s --input <samfile> --design <designfile> --output <output file>

Options:
    -i/--input   [file]  the input SAM/BAM file, default None
    -d/--design  [file]  the input design file, default stdin
    -o/--output  [file]  the output output file, default stdout
    -q/--quality [int]   the minimum quality, default 0
    -h/--help   []      print the help message

""" % {
    "tool": sys.argv[0],
    "commandline": " ".join(sys.argv),
    "message": message
    }
        )
        if error is not None:
            sys.exit(int(error))

    def main():
        """The main program loop."""
        samin = None
        desstream = sys.stdin
        outstream = sys.stdout
        quality = 0
        try:
            shortopt = "i:d:o:q:h"
            longopt = ["input=", "design=", "output=", "quality=", "help"]
            opts, _ = getopt.getopt(sys.argv[1:], shortopt, longopt)
            for opt, ans in opts:
                if opt in ("-i", "--input"):
                    _, ext = os.path.splitext(ans)
                    mode = "rb" if ext == ".bam" else "r"
                    samin = pysam.AlignmentFile(ans, mode)
                elif opt in ("-d", "--design"):
                    desstream = open(ans, "r")
                elif opt in ("-o", "--output"):
                    outstream = open(ans, "w")
                elif opt in ("-q", "--quality"):
                    quality = int(ans)
                elif opt in ("-h", "--help"):
                    usage("Help was asked", error=0)
        except getopt.GetoptError:
            usage("An error occured while parsing commandline parameters", error=100)

        # count the reads per block
        nimbus_count(
            samin,
            bedfile=desstream,
            outstream=outstream,
            quality=quality
        )

        # close the input and output files
        samin.close()
        if desstream is not sys.stdin:
            desstream.close()
        if outstream is not sys.stdout:
            outstream.close()

    #
    main()
