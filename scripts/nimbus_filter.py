#!/bin/env python

import sys
import os
import os.path
import getopt
import pysam

class MismatchCalculator(object):

    def __init__(self):
        pass

    def bases_in_indels(self, cigar=None):
        """
        Determines the number of bases in Indels

        Args:
            cigar:  a list with cigar operations

        Returns:
            2 lists with the with the insertion, deletion and softclipped bases
            (number, number of bases)
        """
        count = [0, 0, 0]
        bases = [0, 0, 0]
        for cig in cigar:
            if cig[0] == 1: #cig[0] == 'I' or
                count[0] += 1
                bases[0] += cig[1]
            elif cig[0] == 2: #cig[0] == 'D' or
                count[1] += 1
                bases[1] += cig[1]
            elif cig[0] == 4: #cig[0] == 'S' or
                count[2] += 1
                bases[2] += cig[1]
        return (count, bases)

    def mismatches(self, nmscore=0, cigar=None):
        """
        Get the number of mismatches based on the NM
        score (Levenstein distance) and the alignment CIGAR.

        Args:
            nmscore: the NM score
            cigar: the CIGAR

        Returns:
            the number of mismatches and the parsed insertion
            and deletion information
        """
        count, bases = self.bases_in_indels(cigar)
        return (
            nmscore - sum(bases),
            count,
            bases
            )


class FilterMismatches(object):

    def __init__(self):
        self.mmc = MismatchCalculator()

    def apply(self, alignment):
        # parse the nm score
        nmscore = None
        tags = alignment.get_tags()
        for tag in tags:
            if tag[0] == "NM":
                nmscore = tag[1]
        #
        if nmscore is None:
            return None

        # get the number of mismatches
        mismatches, count, bases = self.mmc.mismatches(nmscore, alignment.cigartuples)

        # get the total number of discrete differences, where insertions
        #  and deletions are counted as a single difference
        totalscore = mismatches + count[0] + count[1]

        return {
            "total": totalscore,
            "mismatches": mismatches,
            "insertions": count[0],
            "deletions": count[1]
        }


def filter_on_mismatches(samin, samout, discarded=None, maxtotal=6, list_events=False):
    """

    """
    fmms = FilterMismatches()
    for alignment in samin.fetch(until_eof=True):
        if alignment.is_unmapped:
            if list_events:
                sys.stdout.write("%s\tNA\tNA\tNA\tNA\n" % alignment.query_name)
            if discarded is not None:
                discarded.write(alignment)
        else:
            scores = fmms.apply(alignment)
            if scores is not None and scores["total"] <= maxtotal and scores["total"] >= 0:
                samout.write(alignment)
            elif discarded is not None:
                discarded.write(alignment)
            if list_events:
                sys.stdout.write("%s\t%d\t%d\t%d\t%d\n" % (
                    alignment.query_name,
                    scores["total"],
                    scores["mismatches"],
                    scores["insertions"],
                    scores["deletions"]
                    ))
    # end of filter_on_mismatches


if __name__ == "__main__":

    def usage(message="", error=None):
        """Print the usage information for the tool."""
        sys.stdout.write(
"""
%(message)s

Called as: %(commandline)s

Usage:
    %(tool)s --input <inputfile> --output <output file> --paired

Options:
    -i/--input      [file]  the input file, default -
    -o/--output     [file]  the output file, default -
    -d/--discarded  [file]  the output file, default None
    -s/--score      [int]   the maximum total mismatch score of an alignment
    -l/--list       []      list the number of mismatches per alignment in the standard out
    -h/--help   []      print the help message
""" % {
    "tool": sys.argv[0],
    "commandline": " ".join(sys.argv),
    "message": message
    }
        )
        # -p/--paired []      process the reads as pairs (SAM file should be sorted on queryname)
        if error is not None:
            sys.exit(int(error))

    def main():
        # Parse the parameters provided on the commandline
        #
        fnin = "-"
        fnout = "-"
        fndisc = None
        listmm = False
        maxtotalscore = 6
        try:
            # Use getopt to parse the command line
            #
            shortopt = "i:o:hs:d:l"
            longopt = ["input=", "output=", "help", "score=", "discarded=", "list"]
            opts, args = getopt.getopt(sys.argv[1:], shortopt, longopt)
            for opt, ans in opts:
                if opt in ("-i", "--input"):
                    if ans != "-":
                        fnin = os.path.abspath(ans)
                elif opt in ("-o", "--output"):
                    if ans != "-":
                        fnout = os.path.abspath(ans)
                elif opt in ("-d", "--discarded"):
                    fndisc = os.path.abspath(ans)
                elif opt in ("-s", "--score"):
                    maxtotalscore = int(ans)
                elif opt in ("-l", "--list"):
                    listmm = True
                elif opt in ("-h", "--help"):
                    usage("Help was asked", error=0)

        except getopt.GetoptError:
            usage("An error occured while parsing commandline parameters", error=100)

        # Open the alignment in and output files
        #
        bn, ext = os.path.splitext(fnin)
        mode = "rb" if ext == ".bam" else "r"
        infile = pysam.AlignmentFile(fnin, mode)

        # check that the file is in the correct sort order
        #try:
        #    if infile.header['HD']['SO'] != "queryname" and paired:
        #        usage("SAM file is not sorted on queryname", error=101)
        #except KeyError:
        #    pass

        bin, ext = os.path.splitext(fnout)
        mode = "wb" if ext == ".bam" else "w"
        outfile = pysam.AlignmentFile(fnout, mode, template=infile, header=infile.header)

        discfile = None
        if fndisc is not None:
            bin, ext = os.path.splitext(fndisc)
            mode = "wb" if ext == ".bam" else "w"
            discfile = pysam.AlignmentFile(fndisc, mode, template=infile, header=infile.header)

        # Process the data
        #
        filter_on_mismatches(
            samin=infile,
            samout=outfile,
            discarded=discfile,
            maxtotal=maxtotalscore,
            list_events=listmm
        )

        # Close the in and output files after the analysis
        #
        if not infile.closed:
            infile.close()
        if not outfile.closed:
            outfile.close()
        if discfile is not None and not discfile.closed:
            discfile.close()

    main()
