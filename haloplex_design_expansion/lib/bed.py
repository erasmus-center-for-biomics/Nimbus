__author__ = '135001'

import sys

class Bed(object):
    '''
    A class to represent a bed region
    '''
    def __init__(self, chr=None, start=-1, end=-1, strand='*', name='' ):
        self.chr = chr
        self.start = start
        self.end = end
        self.strand = strand
        self.name = name

        if start < 0:
            raise ValueError("start is smaller than 0")
        if end < 0:
            raise ValueError("end is smaller than 0")
        if strand not in ("+", "-", "*"):
            raise ValueError("strand is an unknown value")

    def __eq__( self, other ):
        rval = False
        if other is not None:
            if self.chr == other.chr and self.start == other.start and self.end == other.end and self.strand == other.strand : rval = True
        return rval

    def __str__(self):
        return "%s\t%d\t%d\t%s\t0\t%s" %( self.chr, self.start, self.end, self.name, self.strand)

    def cmp_str(self):
        return "%s\t%d\t%d\t%s" %( self.chr, self.start, self.end, self.strand)

    def __hash__(self):
        return hash( self.cmp_str() )


def bed_from_line( l="" ):
    '''
    Create a BED region from a line
    '''
    r = None
    l = l.rstrip()
    a = l.split( "\t" )
    if len(a) >= 6:
        r = Bed( a[0], int(a[1]), int(a[2]), a[5], a[3] )
    return r

def read_bed( f ):
    '''
    Read the contents of a BED file
    '''
    rval = []

    for l in f:
        b = bed_from_line( l )
        if b is not None:
            rval.append(b)
    # return the
    return rval

def index_bed_file( f=None ):
    """
    Creates a simple index of a bed file

    :param f: an opened file object, supporting the tell method
    :return: a dict with the chromosome as the key and the position of the first occurence
    """
    pchr = None
    lpos = f.tell()
    rval = {}

    # we can not use for l in f because it buffers
    l = f.readline()
    while l != '':

        # move over any headers
        l = l.rstrip()
        a = l.split('\t')
        while len(a) < 3:
            lpos = f.tell()
            l = f.readline()
            a = l.split('\t')

        chr = a[0]
        # we got to the next chromosome, so the
        # end of the previous line was our entry
        if chr != pchr:
            if rval.has_key(chr) : raise ValueError( "%s was already added, so the bed file does not seem to sorted" % chr )
            rval[chr] = lpos

        # set the data for the next cycle
        pchr = chr
        lpos = f.tell()
        l = f.readline()

    # return the index
    return rval

def bediterator( f, index=None, chr=None):
    '''
    Iterates over the entries in a bed file
    :param f: the bed file
    :return: an object of class Bed
    '''

    # indexing code
    if chr is not None:
        if index is None:
            raise ValueError( "index not provided"  )
        if index.has_key( chr ):
            f.seek( index[chr], 0 )
        else:
            raise ValueError( "chromosome '%s' not found in index" % chr )

    # foreach line from the bed file
    # cnt = 0
    l   = f.readline()
    while l != '':
        b = bed_from_line( l )
        # if cnt < 3:
        #     sys.stderr.write( str(l) )
        #     cnt += 1
        if chr is not None and b is not None and b.chr != chr:
            break
        if b is not None:
            yield b
        l = f.readline()