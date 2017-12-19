__author__ = '135001'


import re

'''

Functions to deal with DNA restriction

'''

IUPAC_bases = ('A','C','T','G','R','Y','S','W','K','M','B','D','H','V','N' )

class RestrictionEnzyme(object):
    '''
    A class for
    '''
    def __init__(self, name='', recognition='' ):

        # set the non interpreted data
        self.name = name
        self.recognition_site = recognition

        # make sure the enzyme has a recognition sequence
        assert len(recognition) > 0

        # determine the cut sites for the lagging and lead strand
        self.rcutlead = recognition.replace('|','').find( '/' )
        self.rcutlag  = recognition.replace('/','').find( '|' )
        if self.rcutlag == -1:
            self.rcutlag = self.rcutlead

        # remove the | and signs from the sequence
        tmp = ''.join( [ x for x in self.recognition_site if x in IUPAC_bases ] )
        self.regexp = re.compile( IUPAC2regexp(tmp), flags=re.IGNORECASE )

    def invert(self):
        '''
        Inverts the Restiction enzyme: handy for distance cutter
        '''

        # create the inverse of the recognition sequence
        tmp = []
        for x in self.recognition_site[::-1]:
            if x in IUPAC_bases:
                tmp.append( complement(x) )
            elif x == '/':
                 tmp.append( '|' ) # += '|'
            elif x == '|':
                 tmp.append( '/' ) # += '/'

        # we had only a single cut site
        # if tmp.count('/') == 0:
        #     tmp.replace('|','/')

        #
        nseq = ''.join( tmp )
        if nseq.count('/') == 0:
            nseq = nseq.replace('|','/')

        # get the other
        oth = RestrictionEnzyme( name='inv_' + self.name, recognition=nseq )
        return oth

    def __str__( self ):
        '''

        '''
        return "%s\t%s\t%d\t%d" % (self.name, self.recognition_site, self.rcutlead, self.rcutlag)

def complement(x):
    assert x in IUPAC_bases
    if x == 'A':
        return 'T'
    elif x == 'T':
        return 'A'
    elif x == 'G':
        return 'C'
    elif x == 'C':
        return 'G'
    elif x == 'N':
        return 'N'
    elif x == 'R':
        return 'Y'
    elif x == 'Y':
        return 'R'
    elif x == 'W':
        return 'W'
    elif x == 'S':
        return 'S'
    elif x == 'K':
        return 'M'
    elif x == 'M':
        return 'K'
    elif x == 'B':  # CGT
        return 'V'  # GCT
    elif x == 'D':  # AGT
        return 'H'  # TCA
    elif x == 'H':  # ACT
        return 'D'  # TGA
    elif x == 'V':  # ACG
        return 'B'  # TGC


def palindrome( seq ):
    '''
    Determine whether the enzyme recognizes a palindrome
    (and thus if we can search the genome with a single recognition sequence)

    '''
    tmp = seq.replace('/','')
    tmp = tmp.replace('|','')

    # determine the reverse complement
    cmp = ''.join( [complement(x) for x in tmp[::-1] ] )
    if cmp == tmp :
        return True
    else :
        return False


def IUPAC2regexp( seq ):
    '''
    Decodes a sequence which is IUPAC and convert
    this to a regular expression friendly sequence
    '''
    # convert IUPAC bases
    seq = seq.replace( 'R', '[A,G]' )
    seq = seq.replace( 'Y', '[C,T]' )
    seq = seq.replace( 'S', '[G,C]' )
    seq = seq.replace( 'W', '[A,T]' )
    seq = seq.replace( 'K', '[G,T]' )
    seq = seq.replace( 'M', '[A,C]' )
    seq = seq.replace( 'B', '[C,G,T]' )
    seq = seq.replace( 'D', '[A,G,T]' )
    seq = seq.replace( 'H', '[A,C,T]' )
    seq = seq.replace( 'V', '[A,C,G]' )
    # seq = seq.replace( 'N', '[A,C,T,G]' )
    seq = seq.replace( 'N', '.' )
    # returns the sequence
    return seq

def occurences( sequence, enzyme=None ):
    '''
    Iterates over the recognition sites of the enzyme
    '''
    start = 0

    while True:
        # we are using the buffer interface to avoid copying
        # the sequence string in memory.
        # The finditer solution below does not yield overlapping hits (really annoying!)
        grp = enzyme.regexp.search( buffer(sequence, start) )
        if grp is None :
            break
        else:
            # make sure to increment the search position with the match start + 1
            # without the + 1 we will end up in an infinite loop where the same hit
            # is found over and over again.
            ostart = start
            start += grp.start() + 1

            # yield the start to end of the match
            yield grp.start() + ostart, grp.end() + ostart



if __name__ == '__main__':

    # test a normal cutter
    a = RestrictionEnzyme(name='AluI', recognition='AG/CT')
    print str( a )
    print str( a.invert() )

    # test a distance cutter
    b = RestrictionEnzyme(name='BccI', recognition='CCATCNNNN/N|')
    print str( b )
    print str( b.invert() )

    # test a "weird cutter"
    c = RestrictionEnzyme(name='SfcI', recognition='C/TRYA|G')
    print str( c )
    print str( c.invert() )

    c = RestrictionEnzyme(name='MluCI', recognition='/AATT|')
    print str( c )
    print str( c.invert() )
