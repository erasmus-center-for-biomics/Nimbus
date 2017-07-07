#!/bin/env python

import sys
import math
import operator
import naivetrack.models


__author__ = "R.W.W.Brouwer"
__version__ = "1"
__copyright__ = "Copyright 2015, Erasmus MC"
__credits__ = []
__license__ = "To be decided"
__maintainer__ = "R.W.W.Brouwer"
__email__ = "r.w.w.brouwer@erasmusmc.nl"
__status__ = "testing"


def variant_reader(handle=sys.stdin, chromosome_list=None):
    """
    Generates variants from an input stream
    :param handle:
    :param chromosome_list:
    :yield:
    """
    for line in handle:
        line = line.rstrip()
        variant = naivetrack.models.TrackEntry()
        headers = []
        if len(line) > 0 and line[0] != "#":
            variant.position = naivetrack.models.Position(chromosome_list=chromosome_list)
            variant.position.from_string(line)
            variant.statistics.from_string(line)
            for i in xrange(variant.statistics.number_of_alleles):
                allele = naivetrack.models.Allele()
                allele.from_string(handle.next().rstrip())
                variant.alleles.append(allele)
        else:
            headers.append(line)

        #
        if variant.position is not None:
            yield variant, headers


class VariantParser(object):
    """
    Parses variants from a file stream
    """
    def __init__(self, chromosome_list=None, handle=sys.stdin):
        """
        Initializes the variant parser object and
        initializes the chromosome_list to be shared later

        :param chromosome_list:
        :param handle:
        :return:
        """
        self.chromosome_list = chromosome_list if chromosome_list is not None else []
        self.handle = handle
        self.header = []
        self.reader = variant_reader(handle=self.handle, chromosome_list=self.chromosome_list)

    def __del__(self):
        try:
            if not self.handle.closed:
                self.handle.close()
        except AttributeError:
            pass

    def next(self):
        """
        Get the next variant
        :return: a new variant or None
        """
        v, h = self.reader.next()
        return v, h

    def __iter__(self):
        return self


class BufferedVariantParser(VariantParser):
    """
    A buffered version of the variant parser that keeps the
      last element available
    """

    def __init__(self, **kwargs):
        super(BufferedVariantParser, self).__init__(**kwargs)
        self.last = None

    def next(self):
        v, h = super(BufferedVariantParser, self).next()
        self.last = (v, h)
        return v, h

    def __iter__(self):
        return self


class TrackParser(object):
    """
    Parse variants from multiple files
    """
    def __init__(self, chromosome_list=None):
        self.chromosome_list = chromosome_list if chromosome_list is not None else []
        self.handles = []
        self.bufferedreaders = []
        self.buffer = []
        self.chunksize = 10000
        self.chunki = 0
        self.chunkmax = 500000
        self.curchr = 0
        self.curi = 0

    def __del__(self):
        for h in self.handles:
            if not h.closed:
                h.close()

    def set_maximum_chromosome_size(self, maxsize=None):
        if maxsize is not None:
            self.chunkmax = math.ceil(
                operator.truediv(maxsize, self.chunksize))

    def register_handle(self, handle=sys.stdin):
        """
        register a file handle
        :param handle:
        :return:
        """
        self.handles.append(handle)
        self.bufferedreaders.append(
            [
                BufferedVariantParser(handle=handle, chromosome_list=self.chromosome_list),
                True
            ])

    def next(self):
        """
        Get the next read

        :return: the next variant
        """
        def next_chunk():
            """
            read the next chunk from the input handles
            :return: the maximum and minimum coordinates read
            """
            mincrd = self.chunki * self.chunksize
            maxcrd = (self.chunki + 1) * self.chunksize
            # print "%d - %d" % (mincrd, maxcrd)

            # we need to raise a stopiteration if the streams are
            # processed. However we cannot do this when the last
            # call provides results. So we need to test this
            # before we process results
            nclosed = 0
            for reader in self.bufferedreaders:
                if not reader[1]:
                    nclosed += 1
            if len(self.bufferedreaders) == nclosed:
                raise StopIteration

            # Put data into our temporary buffer
            #
            b = []
            for reader in self.bufferedreaders:
                bufferedreader = reader[0]
                while reader[1]:
                    try:
                        if bufferedreader.last is not None:
                            if bufferedreader.last[0].position.chromosome > self.curchr:
                                break
                            elif bufferedreader.last[0].position.coordinate >= maxcrd:
                                break
                            elif bufferedreader.last[0].position.coordinate < mincrd:
                                break
                            b.append(bufferedreader.last[0])
                        bufferedreader.next()
                    except StopIteration:
                        reader[1] = False
                        break

            b.sort(key=operator.attrgetter("position"))

            # Merge subsequent entries in buffer
            self.buffer = []
            for i in xrange(len(b)):
                if len(self.buffer) > 0:
                    if b[i].position != self.buffer[-1].position:
                        self.buffer.append(b[i])
                    else:
                        self.buffer[-1].add(b[i])
                else:
                    self.buffer.append(b[i])

        # we traversed through our buffer
        #   the while loop lets us traverse
        #   regions without variants
        while self.curi >= len(self.buffer):

            # increase the chunknumber
            self.chunki += 1

            # next chromosome
            if self.chunki > self.chunkmax:
                self.chunki = 0
                self.curchr += 1

            # get the next chunk
            next_chunk()
            self.curi = 0

        rval = self.buffer[self.curi]
        self.curi += 1
        return rval

    def __iter__(self):
        return self
