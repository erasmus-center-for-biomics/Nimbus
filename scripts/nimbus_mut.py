#!/bin/env python

import os
import os.path
import sys
import re
import getopt
import logging
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "../"))
import table
import filters


__version__ = "0.5.dev"


class SampleParser(object):

    def __init__(self):
        self.regexp = re.compile(r"[0-9]*_?(.*)_\[(.*)\]$")

    def get(self, value):
        m = self.regexp.search(value)
        if m is None:
            return None
        else:
            return m.group(1), m.group(2)


class MultiAnnoParser(table.MultiAnno):

    def __init__(self, na_value="-", *args, **kwargs):
        super(MultiAnnoParser, self).__init__(*args, **kwargs)
        self.na_value = na_value
        self.columns = []

    def next(self):
        values = super(MultiAnnoParser, self).next()
        sp = SampleParser()
        tmp = [sp.get(x[0]) for x in values]
        baseinfo = {}
        sampleinfo = {}

        # prepare the header
        if len(self.columns) == 0:
            for i in xrange(len(tmp)):
                if tmp[i] is None:
                    self.columns.append(values[i][0])
                else:
                    if tmp[i][0] not in self.columns:
                        self.columns.append(tmp[i][0])

        # parse the sample information
        for i in xrange(len(tmp)):
            v = tmp[i]

            # set NA values to None
            val = values[i][1]
            if val == self.na_value:
                val = None

            # if a sample column
            if v is not None:
                if v[1] not in sampleinfo.keys():
                    sampleinfo[v[1]] = {}
                sampleinfo[v[1]][v[0]] = val
            else:
                baseinfo[values[i][0]] = val

        # returns the base and sample info
        return baseinfo, sampleinfo


def impute_type(baseinfo):
    """
    impute the variant type from the basic variant information 
    """
    retval = "-"
    if baseinfo['ExonicFunc.refGene'] is not None:
        retval = baseinfo['ExonicFunc.refGene']
    else:
        retval = baseinfo['Func.refGene']
    if baseinfo['Func.refGene'] == "splicing":
        retval = "splicing"
    return retval


def multianno_to_mut(hin, hout, columnfilters=None):
    """
    Convert a multianno table to a mut
    """
    parser = MultiAnnoParser(h=hin)
    writer = table.Writer(h=hout)
    writer.columns = ["chr", "start", "end", "sample", "type"]

    # foreach entry in the input file
    for baseinfo, sampleinfo in parser:
        vartype = impute_type(baseinfo)

        # append columns on the first entry
        if len(writer.columns) == 5:
            writer.columns.extend(parser.columns[3:])

        # add one entry per sample to the linebuffer
        linebuffer = []
        for sample, info in sampleinfo.items():
            data = [baseinfo["Chr"], baseinfo["Start"], baseinfo["End"],
                    sample, vartype]

            # append the column information
            for name in writer.columns[5:]:
                if name in baseinfo.keys():
                    data.append(baseinfo[name])
                elif name in info.keys():
                    data.append(info[name])
                else:
                    data.append(None)
            linebuffer.append(dict(zip(writer.columns, data)))

        # filter entries
        if columnfilters is not None:
            fltbuffer = []
            for data in linebuffer:
                keep = 0
                for cflt in columnfilters:
                    try:
                        cflt.check(data)
                        keep += 1
                    except filters.RejectedException:
                        pass
                    except LookupError:
                        pass
                if keep == len(columnfilters):
                    fltbuffer.append(data)
        else:
            fltbuffer = linebuffer

        # print the entries that are not filtered data
        for data in fltbuffer:
            writer.add(indict=data)

    # log the data considered at the end
    logging.info("Read %d lines from the input handle" % parser._linecount)
    logging.info("Wrote %d lines to the output handle" % writer._linecount)


if __name__ == "__main__":
    def usage(message="", error=None):
        """
        Prints the usage information and messages
        """
        values = {"message": message, "commandline": " ".join(sys.argv),
                  "tool": sys.argv[0], "version": __version__}

        sys.stderr.write("""
%(message)s

Called as: %(commandline)s

Usage:
    %(tool)s --input [input.multianno.txt] --out [output.mut] --filter 'Depth_frequency:<0.2,->'

Options:
    -i/--input  [file]      the input multianno file from Annovar
    -o/--output [file]      the output mut file (see
                            http://software.broadinstitute.org/software/igv/MUT
                            for format specification)
    -f/--filter [string]    a string on which to filter the mut output file.
                            These strings have the following format 
                            COLUMN:<minimum,maximum> for a numeric filter (- for undefined)
                            of COLUMN:VALUE for an equals filter

version: %(version)s
""" % values)
        if error is not None:
            sys.exit(int(error))

    logging.basicConfig(
            level=logging.INFO,
            format='[%(asctime)s %(name)s %(process)d %(levelname)-8s]: %(message)s')

    # initialize the variables
    fin = sys.stdin
    fout = sys.stdout
    factory = filters.FilterFactory()
    allfilters = []

    # parse the commandline options
    shortopt = "i:o:f:h"
    longopt = ["input=", "output=", "filter=", "help"]
    try:
        opts, argus = getopt.getopt(sys.argv[1:], shortopt, longopt)
        for o, a in opts:
            if o in ("-i", "--input"):
                fin = open(os.path.abspath(a), 'r')
            elif o in ("-o", "--output"):
                fout = open(a, 'w')
            elif o in ("-f", "--filter"):
                flt = factory.get_dict_filter(a)
                if flt is None:
                    usage("Could not parse a filter from %s" % a)
                allfilters.append(flt)
            elif o in ('-h', "--help"):
                usage("Help was asked", error=0)
    except getopt.GetoptError as err:
        usage(str(err), error=101)

    # run the main loop
    multianno_to_mut(fin, fout, allfilters)

