#!/bin/env python

import os
import os.path
import sys
import re
import getopt
from operator import itemgetter, truediv

# add the library path relative to the script location
libpath = os.path.abspath(os.path.join( os.path.dirname(sys.argv[0]), '../' ))
libpath = os.path.join(libpath, 'lib')
sys.path.append(libpath)

# load our custom libraries
from bed import Bed, bediterator, index_bed_file
from FastA import fasta_chromosome, index_fasta
from dna_restriction import complement
from wordtree import WordTree


def usage(mess="", error=None):
    sys.stdout.write("""
%s

called as: %s

Usage:
%s -b <bed> -a <all>-f <fasta file> -o <out.bed>

Options:
-b/--bed    A bed file with the amplicons to check
-a/--all    A bed file with all the possible restriction fragments on that chromosome
-o/--out    The output file
-h/--help   prints this message

    """ % (mess, ' '.join(sys.argv), sys.argv[0]))

    if error is not None:
        sys.exit(int(error))


def distance(a, b, max=None):
    rval = len(a)

    for i in xrange(len(a)):
        if a[i] == b[i]:
            rval += 1
        if max is not None and rval > max:
            break
    return rval

def getBedKey(sequence, bed, ksize=7):
    '''
    This function gets a sequence key for the current bed entry

    Let our sequence be

    GATG AAGC TGGA
    (1)        (2)
    and our key size is 4

    ---------------
    * for a forward bed entry the 2 sub keys will be

     (GATG) and rc(TGGA):TCCA
       (1)        (2)

    ---------------
    * for a reverse bed entry the keys are

    rc(TGGA):TCCA and GATG
      (2)           (1)

    Rules:
        if +
            (1) + rc(2)
        else:
            rc(2) + 1

    Evidence:
        the reverse complement of our sequence is

        TCCA GCTT CATC
        (1)'      (2)'
        ---------------
        * the + strand case in the last example is now the - strand case, so lets take the
          rules we defined for the - strand case

        rc(CATC):GATG   equals the first forward key in the original sequence

        and

        TCCA             equal the second forward key in the original sequence

        ---------------
        * the - strand case of the last example is now + strand case, so

        TCCA            equals the first reverse key in the original sequence

        and

        rc(CATC):GATG   does not equal the second reverse key in the original sequence

    '''
    seqa = sequence[bed.start : (bed.start + ksize)]
    seqb = sequence[(bed.end - ksize) : bed.end]
    key  = None
    if bed.strand == '+' :
        tmp = [complement(b) for b in seqb[::-1]]
        key = seqa + '|' + ''.join(tmp)
    else:
        tmp = [complement(b) for b in seqb[::-1]]
        key = ''.join(tmp) + '|' + seqa

    # return the key
    return key


if __name__ == '__main__':

    fnin = None
    fnall = None
    fout = sys.stdout
    fnfasta = None
    ksize = 10
    mismatch = 2
    sbases = 5
    verbose = False
    try:
        opts, args = getopt.getopt(sys.argv[1:], 'b:a:o:f:k:m:s:hv',  ['all=', 'bed=', 'out=', 'fasta=', 'ksize=', 'mismatches=', 'safebases=', 'help', 'verbose'])

        for o,a in opts:
            if o in ('-b','--bed'):
                fnin = os.path.abspath(a)
            elif o in ('-a','--all'):
                fnall = os.path.abspath(a)
            elif o in ('-v','--verbose'):
                verbose = True
            elif o in ('-f', '--fasta'):
                fnfasta = os.path.abspath(a)
            elif o in ('-k', '--ksize'):
                ksize = int(a)
            elif o in ('-m', '--mismatches'):
                mismatch = int(a)
            elif o in ('-s', '--safebases'):
                sbases = int(a)
            elif o in ('-o','--out'):
                fout = open(os.path.abspath(a), 'w')
            elif o in ('-h', '--help'):
                usage( 'Help was asked', error=1 )
    except getopt.GetoptError, err:
        usage(str(err), error=2)

    #
    if fnfasta is None or not os.path.exists(fnfasta):
        usage( "FastA file %s does not exist" % fnfasta, error=3 )

    # opens the fasta file
    ffasta = open(fnfasta, 'rU')
    sys.stderr.write( "Indexing FastA file %s\n" % fnfasta )
    fai = index_fasta(ffasta)

    #
    # We will keep the keys for the amplicon design in memory and
    # match the entries in the complete restriction map against them
    #

    # generate a set with all the keys we are searching for
    tree = WordTree()
    fin = open(fnin, 'rb')
    bedidx = index_bed_file( fin )      # input file should be sorted on chromosome
    sys.stderr.write( "Determining all keys in the design\n" )

    # for each chromosome in fasta
    cnt = 0
    for chr in sorted(fai.keys()):

        if not chr in bedidx.keys() or not chr in fai.keys():
            continue

        # get all the Bed entries for the current chromosome
        bed  = [b for b in bediterator(fin, index=bedidx, chr=chr)]

        # get the sequence and sequence id for the current chromosome
        seqid, sequence = fasta_chromosome(ffasta, index=fai, chr=chr )
        sequence = sequence.upper()

        # get all the key sequences
        for b in bed:
            k = getBedKey( sequence, b, ksize=ksize)
            if verbose:
                sys.stderr.write( "%s\t%s\n" % (str(b), k) )
            tree.add_word( k )
            cnt += 1

    sys.stderr.write("Added %d words to the Tree (%d)\n" % (tree._words_added, cnt))

    # close the input file
    fin.close()

    # process the design
    sys.stderr.write("Aggregating entries\n")
    fall = open(fnall, 'rb')
    bedidx = index_bed_file(fall)

    processed = 0

    # set the protected zone and maximum number of mismatches per key
    tree._protected = sbases
    tree._maxscore = mismatch
    tree._reset_char = '|'

    # for each chromosome in fasta file
    for chr in sorted(fai.keys()):

        if not chr in bedidx.keys() or not chr in fai.keys():
            continue

        # state when we are switching chromosomes
        sys.stderr.write("\tProcessing sequence %s for %s\n" % (chr,fnall))

        # get the sequence and sequence id for the current chromosome
        seqid, sequence = fasta_chromosome(ffasta, index=fai, chr=chr)
        sequence = sequence.upper()

        # get all the Bed entries for the current chromosome
        for b in bediterator(fall, index=bedidx, chr=chr):

            # check the presence of the key in design keyset
            key  = getBedKey(sequence, b, ksize=ksize)

            # tell us how far we are
            if processed % 1000000 == 0:
                sys.stderr.write("\t\t%d regions processed\n" % processed)

            # if the key is present in the tree (with a maximum number of $mismatch mismatches)
            p = tree.is_present(key)
            if p == True:
                b.name += ';' + key
                fout.write("%s\n" % str(b))

            # increase the process counter
            processed += 1

    # close the input and output files
    fall.close()
    ffasta.close()
    fout.close()

