#!/bin/env python

import sys
import argparse
import gzip


if __name__ == "__main__":

    def main():
        """The main program loop."""
        parser = argparse.ArgumentParser(
            prog=sys.argv[0],
            description="""
            A script to merge the variant calls from 
            multiple chromosomes to a single var file.
            """)
        parser.add_argument(
            "--output",
            dest="output",
            help="The output file",
            type=str, nargs="?", default="merged.var")
        parser.add_argument(
            "--input",
            dest="input",
            help="The input files",
            type=str, nargs="+", default=[])
        args = parser.parse_args()

        # set the output stream
        outstream = sys.stdout
        if args.output.endswith(".gz"):
            outstream = gzip.open(args.output, "wb")
        else:
            outstream = open(args.output, "w")

        # foreach input file
        for fname in args.input:

            # open it and copy all data lines
            fin = open(fname, "r") if not fname.endswith(".gz") else gzip.open(fname, "rb")
            for line in fin:
                if line.startswith("#"):
                    continue
                outstream.write(line)
            fin.close()

        # close the output file
        if outstream is not sys.stdout:
            outstream.close()
    #
    main()
