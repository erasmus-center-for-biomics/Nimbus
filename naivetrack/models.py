#!/bin/env python

# library imports
import re

__author__ = "R.W.W.Brouwer"
__version__ = "1"
__copyright__ = "Copyright 2015, Erasmus MC"
__credits__ = []
__license__ = "To be decided"
__maintainer__ = "R.W.W.Brouwer"
__email__ = "r.w.w.brouwer@erasmusmc.nl"
__status__ = "testing"

# The pattern with which to split alleles.
#
#  We will be doing this a lot, so we have
#   compiled the pattern beforehand
allele_tag_pattern = re.compile(", ?")


# 
AlleleFactors = {
    "Allele": 0,
    "Sample": 1,
    "Strand": 2
}


class Allele(object):
    """
    A class to represent an allele.

    Each allele is characterized by a sequence with
    its cumulative depth and quality. An allele has an
    info field used to hold information as tags (like sample).
    """

    # set the slots of the object
    __slots__ = ['depth', 'quality', 'sequence', 'tags']

    def __init__(self, depth=0, quality=0, sequence='', tags=None):
        self.depth = depth
        self.quality = quality
        self.sequence = sequence.upper()
        self.tags = tags if tags is not None else []

    def from_string(self, s=""):
        """
        Loads an allele from a string
        :param s: the string to process
        :return:
        """
        a = s.split('\t')
        if len(a) != 5 or a[0] != "":
            raise ValueError("Could not retrieve allele from string '%s'" % s)
        self.sequence = a[1]
        self.depth = int(a[2])
        self.quality = int(a[3])
        self.tags = allele_tag_pattern.split(a[4])

    def attribute(self, value=0):
        """
        Gets an attribute from the allele. The attributes are
          the concatenation of the sequence and the tags.

        :param value:
        :return:
        """
        assert value >= 0
        assert value <= len(self.tags)

        if value == 0:
            return self.sequence
        else:
            return self.tags[value-1]

    def add_tag(self, value=""):
        """
        Adds the value @value to the current allele

        :param value: the value to be added
        :return: None
        """
        self.tags.append(value)

    def __str__(self):
        return "\t%(sequence)s\t%(depth)d\t%(quality)d\t%(tags)s" % {
            "sequence": self.sequence, "depth": self.depth, "quality": self.quality, "tags": ", ".join(self.tags)}

    def __cmp__(self, other):
        """
        Compares the current allele to another.

        First the alleles will be compared on sequence. If this is equal, the tags are compared.

        :param other: the other allele
        :return: -1, 0, 1
        """
        if self.sequence < other.sequence:
            return -1
        elif self.sequence > other.sequence:
            return 1
        else:
            if len(self.tags) < len(other.tags):
                return -1
            elif len(self.tags) > len(other.tags):
                return 1
            else:
                for i in xrange(len(self.tags)):
                    if self.tags[i] < other.tags[i]:
                        return -1
                    elif self.tags[i] > other.tags[i]:
                        return 1
        # if we got here, both are equal
        return 0


class Position(object):
    """
    A class to represent a genome position

    Genome positions have a chromosome, a coordinate and a reference base. Reference bases
    in genome positions can be upper or lowercase. With the from_string method, Positions can
    be initialized from a string.
    """

    # set the slots of the object
    __slots__ = ['chromosome', 'coordinate', 'sequence', '__chromosome_list']

    def __init__(self, chromosome="", coordinate=-1, sequence="", chromosome_list=None):
        """
        Initializes a new object

        :param chromosome: the chromosome
        :param coordinate: the coordinate
        :param sequence: the sequence
        :return: the constructed object
        """
        self.chromosome = -1
        self.__chromosome_list = chromosome_list
        self.__set_chromosome(chromosome)
        self.coordinate = coordinate
        self.sequence = sequence

    def __set_chromosome(self, chromosome=""):
        if self.__chromosome_list is not None:
            try:
                self.chromosome = self.__chromosome_list.index(chromosome)
            except ValueError:
                self.__chromosome_list.append(chromosome)
                self.chromosome = self.__chromosome_list.index(chromosome)
        else:
            self.__chromosome_list = [chromosome]
            self.chromosome = 0

    def from_string(self, s=""):
        """
        Initializes a position from a string
        :param s: a string representing a position
        :return:
        """
        a = s.split(':')
        if len(a) >= 3:
            self.__set_chromosome(a[0])
            self.coordinate = int(a[1])
            self.sequence = a[2]
        else:
            a = s.split('\t')
            if len(a) >= 3:
                self.__set_chromosome(a[0])
                self.coordinate = int(a[1])
                self.sequence = a[2]
            else:
                raise ValueError("Cannot obtain Position object from string '%s'" % s)

    def get_chromosome(self):
        """
        Get the chromosome
        :return: the chromosome as a string
        """
        return self.__chromosome_list[self.chromosome]

    def __str__(self):
        """
        Converts the current entry to a string
        :return:
        """
        return "%(chromosome)s\t%(coordinate)d\t%(sequence)s" % {
            "chromosome": self.get_chromosome(), "coordinate": self.coordinate, "sequence": self.sequence}

    def __cmp__(self, other):
        """
        Compares 2 positions

        :param other: the other object
        :return: -1,0,1
        """
        if self.chromosome < other.chromosome:
            return -1
        elif self.chromosome > other.chromosome:
            return 1
        else:
            if self.coordinate < other.coordinate:
                return -1
            elif self.coordinate > other.coordinate:
                return 1
            else:
                return 0


class Statistics(object):
    """
    An object to hold the position statistics at a single position of the genome
    """
    __slots__ = ['number_of_alleles', 'depth', 'quality']

    def __init__(self, numberofalleles=0, depth=0, quality=0):
        """
        Initializes the statistics object
        :param numberofalleles:
        :param depth:
        :param quality:
        :return:
        """
        self.number_of_alleles = numberofalleles
        self.depth = depth
        self.quality = quality

    def from_string(self, s=""):
        """
        Parse the summary statistics from a string

        :param s: the string to parse the statistics from
        :return:
        """
        a = s.split('\t')
        if len(a) >= 6:
            self.number_of_alleles = int(a[3])
            self.depth = int(a[4])
            self.quality = int(a[5])

    def add_statistics(self, other):
        """
        Adds another statistics to the current statistic
        :param other:
        :return:
        """
        self.number_of_alleles += other.number_of_alleles
        self.depth = self.depth
        self.quality = self.quality

    def add_allele(self, allele):
        """
        Adds a allele to the current statistics
        :param allele: the allele to add
        :return: none
        """
        self.number_of_alleles += 1
        self.depth += allele.depth
        self.quality += allele.quality

    def __str__(self):
        return "%(number_of_alleles)d\t%(depth)d\t%(quality)d" % {
            "number_of_alleles": self.number_of_alleles, "depth": self.depth, "quality": self.quality}


class TrackEntry(object):
    """
    A class to represent an entry in the variant file format
    """
    __slots__ = ['position', 'alleles', 'statistics']

    def __init__(self, position=None, alleles=None, statistics=None):
        """
        Initialize the variant
        :param position:
        :param alleles:
        :param statistics:
        :return:
        """
        self.position = position
        self.alleles = alleles if alleles is not None else []
        self.statistics = statistics if statistics is not None else Statistics()

    def __str__(self):
        return "%(position)s\t%(statistics)s\n%(alleles)s" % {
            "position": str(self.position),
            "statistics": str(self.statistics),
            "alleles": '\n'.join([str(a) for a in self.alleles])
        }

    def recalculate_statistics(self):
        """
        Recalculate the statistics for the currently assigned alleles
        :return:
        """
        self.statistics = Statistics()
        for a in self.alleles:
            self.statistics.add_allele(a)

    def add(self, other):
        """
        Add another variant to the current
        :param other:
        :return:
        """
        for allele in other.alleles:
            self.add_allele(allele)

    def add_allele(self, allele=None):
        """
        Adds an allele to the current variants
        :param allele: the allele to add
        :return:
        """
        self.alleles.append(allele)
        self.statistics.add_allele(allele)

    def filter_alleles(self, keep=None):
        """
        Filters the alleles in-place
        :param keep: a True/False list with the alleles to keep
        :return:
        """
        assert keep is not None
        assert len(keep) != len(self.alleles)
        self.alleles = [self.alleles[i] for i in xrange(len(keep)) if keep[i]]

    @staticmethod
    def exact_match(x="", y=""):
        return x == y

    @staticmethod
    def contains_match(x="", y=""):
        return x in y

    def allele_sequence(self, seq="", exact=True):
        """
        Mark all the alleles with the sequence seq
        :param seq: the sequence to match
        :param exact: should we do an exact sequence match or a contains?
        :return: a list with True/False with the same length as self.alleles. True
        """
        matcher = self.exact_match if exact else self.contains_match
        rval = [matcher(seq, a.sequence) for a in self.alleles]
        return rval

    def is_not_reference(self):
        """
        Determine which alleles are not the reference base
        :return: a list with True/False
        """
        isref = self.allele_sequence(self.position.sequence.upper(), exact=True)
        notref = [not x for x in isref]
        return notref

    def n_reference(self):
        """
        Determine the number of reference alleles
        :return: the number of reference alleles
        """
        n = 0
        u = self.position.sequence.upper()
        for a in self.alleles:
            if a.sequence == u:
                n += 1
        return n


if __name__ == "__main__":
    pass
