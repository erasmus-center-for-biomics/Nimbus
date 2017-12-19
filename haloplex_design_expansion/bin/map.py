#!/bin/env python

import os
import os.path
import sys
import re
import getopt
from operator import itemgetter

# add the library path relative to the script location
libpath = os.path.abspath( os.path.join( os.path.dirname(sys.argv[0]), '../' ) )
libpath = os.path.join( libpath, 'lib' )
sys.path.append( libpath )

# load our custom libraries
from FastA import FastAIterator
from dna_restriction import RestrictionEnzyme, occurences, complement, palindrome


def recognition_sequences( sequence, enzymes=[] ):
    '''
    Process the chromosome sequence with the restriction enzymes and
    determine all the sites recognized by the enzymes
    '''
    # get all the locations
    locations = []
    for enzyme in enzymes:
        sys.stderr.write("\tenzyme %s" % enzyme.name)
        cnt = 0
        for oc in occurences(sequence,enzyme):
            locations.append([oc[0], oc[1], enzyme.name])
            cnt += 1
        # report the number of sites
        sys.stderr.write("\t%d sites\n" % cnt)

    # sort the locations on start position
    locations.sort(key=itemgetter(0))
    return locations


def scan_genome(fasta=FastAIterator(), enzymes=[], minwidth=30, maxwidth=600, fout=sys.stdout):
    '''
    Scans the genome for the restriction sites
    '''
    #rval = []

    enzymes_lst = dict( (e.name, e) for e in enzymes )

    # process each chromosome individually
    for id, seq in fasta:
        sys.stderr.write("Scanning sequence %s\n" % id)

        # get all the enzyme recognition sites
        fnd = recognition_sequences(seq, enzymes)
        cnt = 0

        # process the forward strand
        # determine all the sites where the enzyme cut the reference sequence
        cutsites = []
        for f in fnd:
            cutsites.append((f[0] + enzymes_lst[f[2]].rcutlead, f[2]))
        cutsites.sort(key=itemgetter(0))

        prev = None
        for loc in cutsites:
            if prev is not None:
                width = loc[0] - prev[0]

                # make sure we should interpre the
                if width >= minwidth and width <= maxwidth:

                    # determine the + strand restriction fragment
                    entry = (id, prev[0], loc[0], prev[1] + "-" + loc[1], "+")
                    fout.write("%s\t%d\t%d\t%s\t0\t%s\n" % entry)
                    cnt += 1
            prev = loc

        # process the reverse strand
        cutsites = []
        for f in fnd:
            cutsites.append((f[0] + enzymes_lst[f[2]].rcutlag, f[2]))
        cutsites.sort(cutsites, key=itemgetter(0))

        prev = None
        for loc in cutsites:
            if prev is not None:
                width = loc[0] - prev[0]

                # make sure we should interpre the
                if width >= minwidth and width <= maxwidth:

                    # determine the + strand restriction fragment
                    entry = (id, prev[0], loc[0], prev[1] + "-" + loc[1], "-" )
                    fout.write("%s\t%d\t%d\t%s\t0\t%s\n" % entry)
                    cnt += 1
            prev = loc

        #
        sys.stderr.write("Found %d fragments between %d and %d base pairs\n" % (cnt, minwidth, maxwidth))


def load_restriction_enzymes(f=sys.stdin):
    '''
    loads a restriction enzym mix from a file
    '''
    rval = []
    for l in f:
        l = l.rstrip()
        a = l.split('\t')
        if l[0] != '#' and len(a) == 2:
            r = RestrictionEnzyme(a[0], a[1])
            rval.append(r)
            if not palindrome(r.recognition_site):
                rval.append(r.invert())

    # returns the mixes
    return rval


def usage( mess="", error=None ) :
    sys.stdout.write("""
%s

called as: %s

Usage:
%s -e <enzymes> -f <fasta file> -o <out.bed>

Options:
-e/--enzymes    the restriction enzyme mix used
                - format name\tsite
                - the site should be as an IUPAC code with / to
                  indicate the cut site in the leading strand and
                  | for the cut site in the lagging strand. For blunt
                  cutters only indicate the leading strand
-f/--fasta      the fasta genome
-o/--out        the output file
-h/--help   prints this message

    """ % (mess, ' '.join(sys.argv), sys.argv[0]))
    if error is not None:
        sys.exit( int(error) )


if __name__ == '__main__':

    #
    fenzyme = sys.stdin
    fnfasta = None
    fout = sys.stdout

    try:
        opts, args = getopt.getopt(sys.argv[1:], 'e:f:o:h', ['enzymes=', 'fasta=', 'out=', 'help'])
        for o,a in opts:
            if o in ('-e','--enzymes'):
                fenzyme = open( os.path.abspath(a), 'r' )
            elif o in ('-o','--out'):
                fout = open( os.path.abspath(a), 'w' )
            elif o in ('-f', '--fasta'):
                fnfasta = os.path.abspath(a)
            elif o in ('-h', '--help'):
                usage( 'Help was asked', error=1 )

    except getopt.GetoptError, err:
        usage(str(err), error=2)

    #
    if fnfasta is None or not os.path.exists(fnfasta):
        usage("FastA file %s does not exist" % fnfasta, error=3)

    # load the enzymes
    enzymes = load_restriction_enzymes( fenzyme )
    fenzyme.close()
    sys.stderr.write('Found %d enzymes\n' % len(enzymes))

    # prepare the fasta iterator
    ffasta = open( fnfasta, 'rU' )
    fasta  = FastAIterator( f=ffasta )
    sys.stderr.write('Loading %s as reference\n' % fnfasta )

    # get all the recognition sites
    scan_genome( fasta=fasta, enzymes=enzymes, fout=fout )
    ffasta.close()

    # write the output to a BED file
    #for s in sites:

    fout.close()
