__author__ = '135001'


'''

A library to work with FastA files from python

'''


import sys
import os
import os.path

usepysam = False

try:
    import pysam
    usepysam = True
except ImportError as e:
    pass
    # raise ImportError(e)

def index_fasta( f ):
    """
    Makes a simple index for a FastA file
    :param f: an opened fasta file
    :return: a dict with the chromosomes as keys and the tell positions as values
    """
    lpos = f.tell()
    rval  = {}

    # for each line in the fasta
    l  = f.readline()
    while l != "" :
        # remove the endline
        l = l.rstrip()

        # add the previous tell to the dict if we found a new identifier
        if l[0] == '>':
            rval[ l[1:] ] = lpos

        # prepare for the next line
        lpos = f.tell()
        l = f.readline()

    # returns the index
    return rval

def fasta_chromosome( f, chr=None, index=None ):
    """
    gets the sequence for a chromosome
    :param f: the opened fasta file
    :param chr: the chromosome to get
    :param index: the index to point to the start
    :return: a sequence id, sequence tuple
    """
    assert index is not None
    assert chr is not None

    if index.has_key( chr ):
        f.seek( index[chr] )
    else:
        raise ValueError( "chromosome %s not found in index" % chr )

    seqid = f.next().rstrip()
    if seqid[0] != '>':
        raise ValueError( "The index does not seem to correspond to the file" )
    seq   = ''
    for l in f:
        if l[0] != '>':
            seq += l.rstrip()
        else:
            break
    return seqid[1:], seq


def FastAIterator( f=sys.stdin, index=None, chr=None ):
    '''
    Iterate over the sequences in a FastA file
    '''

    id = None
    seq = None

    # process the lines in the output
    for l in f:
        # remove the end line sign
        l = l.rstrip()

        # if we found a new header
        if l[0] == '>':

            # only output a sequence if we got an Id
            if id is not None:
                yield id, seq

            # iterator wait-point here
            id  = l[1:]
            seq = ''

        if len(l) > 0 and l[0] != '>':
            seq += l

    # output the last chromosome
    if id is not None:
        yield id, seq

class FastA(object):

    def __init__(self, fn_fasta=None ):

        if fn_fasta == None:
            raise ValueError( "fn_fasta can not be None" )
        elif not os.path.exists( str(fn_fasta) ) :
            raise ValueError( "FastA file %s was not found" % fn_fasta )

        #
        self._fasta    = fn_fasta
        self._usepysam = usepysam
        if self._usepysam:
            #
            self._pysamfasta = pysam.Fastafile( self._fasta )

    def get_sequence(self, rid=None, start=-1, end=-1, length=-1 ):
        '''

        '''
        assert rid != None
        assert start > 0
        if end == -1: end = start + length
        if length == -1: length = end - start
        assert end > 0
        assert length >= 0

        if self._usepysam:
            return self._pysam_get_sequence( rid=rid, start=start, end=end )
        else:
            raise NotImplementedError( "We have not yet implemented a pysam independent method for extracting sequences from a fasta file" )

    def _pysam_get_sequence( self, rid=None, start=-1, end=-1 ) :
        '''


        '''
        return self._pysamfasta.fetch( reference=rid, start=start, end=end )



if __name__ == '__main__':
    pass
