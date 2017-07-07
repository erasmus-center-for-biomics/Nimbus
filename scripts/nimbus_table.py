import os
import os.path
import sys
import getopt
import logging
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "../"))
import table

__version__ = "1.0"

def multianno_to_table(fin=sys.stdin, fout=sys.stdout):
    """Convert a multianno file to a regular table."""

    # prepare the in and output tables
    parser = table.MultiAnno(h=fin)
    writer = table.Writer(h=fout)
    writer.na_value = "NA"

    # add each row to the new file
    for row in parser:
        if len(writer.columns) == 0:
            writer.columns.extend(parser.header)
        writer.add(tuplelist=row)


if __name__ == "__main__":

    def usage(message="", error=None):
        """Print the usage information and messages."""
        values = {"message": message, "commandline": " ".join(sys.argv),
                  "tool": sys.argv[0], "version": __version__}

        sys.stderr.write("""
%(message)s

Called as: %(commandline)s

Usage:
    %(tool)s --input [input.multianno.txt] --out [output.txt] 

Options:
    -i/--input  [file]      the input multianno file from Annovar
    -o/--output [file]      the output table

version: %(version)s
""" % values)
        if error is not None:
            sys.exit(int(error))

    def main():
        """The main loop."""
        logging.basicConfig(
            level=logging.INFO,
            format='[%(asctime)s %(name)s %(process)d %(levelname)-8s]: %(message)s'
        )

        # initialize the variables
        fin = sys.stdin
        fout = sys.stdout

        shortopt = "i:o:h"
        longopt = ["input=", "output=", "help"]
        try:
            opts, argus = getopt.getopt(sys.argv[1:], shortopt, longopt)
            for o, a in opts:
                if o in ("-i", "--input"):
                    fin = open(os.path.abspath(a), 'r')
                elif o in ("-o", "--output"):
                    fout = open(a, 'w')
                elif o in ('-h', "--help"):
                    usage("Help was asked", error=0)
        except getopt.GetoptError as err:
            usage(str(err), error=101)

        # convert the multianno file to a normal tab-
        #  delimited table
        multianno_to_table(fin, fout)

        # close the in and output files
        if not fin.closed and fin is not sys.stdin:
            fin.close()
        if not fout.closed and fout is not sys.stdout:
            fout.close()
    #
    main()
