#!/bin/env python

# Python standard libraries
import operator
import naivetrack.models
import naivetrack.utils


__author__ = "R.W.W.Brouwer"
__version__ = "1"
__copyright__ = "Copyright 2016, Erasmus MC"
__credits__ = []
__license__ = "To be decided"
__maintainer__ = "R.W.W.Brouwer"
__email__ = "r.w.w.brouwer@erasmusmc.nl"
__status__ = "testing"


class DisqualifiedSampleException(Exception):
    """
    an exception to tell that the sample is finished
    """
    pass


class RangeCheck(object):
    """

    """

    def __init__(self, lower=None, upper=None):
        self.lower = lower
        self.upper = upper

    def __call__(self, value):
        """
        Performs a range check for value.

        :param value: the value to check
        :return: False if the value fails the range check. True
          if the value passes the range check.
        """
        if self.lower is not None:
            if self.__is_lt__(value):
                return False
        if self.upper is not None:
            if self.__is_gt__(value):
                return False
        return True

    def __is_lt__(self, value):
        return True if value <= self.lower else False

    def __is_gt__(self, value):
        return True if value >= self.upper else False

    def __str__(self):
        return "RangeCheck(" + str(self.lower), ", " + str(self.upper) + ")"


class SampleFilter(object):

    def __init__(self):
        """

        :return:
        """
        self.factor_attributes = [
            naivetrack.models.AlleleFactors["Allele"],
            naivetrack.models.AlleleFactors["Sample"]]
        self.sample_attributes = [
            naivetrack.models.AlleleFactors["Sample"]]
        self.value_filters = {}
        self.frequency_filters = {}

        # filter settings
        self.remove_N_alleles = False
        self.remove_all_reference = False

    def n_filters_defined(self):
        n = 0
        n += len(self.value_filters)
        n += len(self.frequency_filters)
        if self.remove_N_alleles:
            n += 1
        if self.remove_all_reference:
            n += 1
        return n

    def __call__(self, variant=None, *args, **kwargs):
        """

        :param variant:
        :param args:
        :param kwargs:
        :return:
        """

        # calculate the reference sequence once
        refseq = variant.position.sequence.upper()

        # actual start of the function
        assert variant is not None
        toremove = []
        outvar = naivetrack.models.TrackEntry(position=variant.position)

        # if we want to remove N
        if self.remove_N_alleles:
            for a in variant.alleles:
                if 'N' in a.sequence and a.sequence not in toremove:
                    toremove.append(a.sequence)

        #
        cache = [None, None]

        # if we have value filters
        if len(self.value_filters) > 0 or len(self.frequency_filters) > 0:
            #
            cache[0] = naivetrack.utils.aggregate_allele_information(variant, attributes=self.factor_attributes)
            if len(self.frequency_filters) > 0:
                cache[1] = naivetrack.utils.aggregate_allele_information(variant, attributes=self.sample_attributes)

            # foreach allele to check
            for altseq, values in cache[0].items():
                if altseq == refseq:
                    continue
                vote = 0
                for sample, data in values.items():
                    try:
                        # check the qualifiers in series
                        for qualifier, check in self.value_filters.items():
                            try:
                                if not check(data[qualifier]):
                                    raise DisqualifiedSampleException()
                            except KeyError:
                                # if we could not get the key, our range check failed
                                raise DisqualifiedSampleException()

                        for qualifier, check in self.frequency_filters.items():
                            try:
                                f = operator.truediv(data[qualifier], cache[1][sample][qualifier])
                                if not check(f):
                                    raise DisqualifiedSampleException()
                            except KeyError:
                                # if we could not get the key, our range check failed
                                raise DisqualifiedSampleException()
                            except ZeroDivisionError:
                                raise DisqualifiedSampleException()
                    except DisqualifiedSampleException:
                        vote += 1
                # if all the samples failed the checks, remove the alleles with this specific alternate sequence
                if vote >= len(values) and altseq != refseq:
                    toremove.append(altseq)

        # remove the marked alleles
        for a in variant.alleles:
            if a.sequence not in toremove:
                outvar.add_allele(a)

        if self.remove_all_reference and outvar.n_reference() == len(outvar.alleles):
            outvar = None
        return outvar

    def add_value_filter(self, qualifier=None, lowerbound=None, upperbound=None):
        """
        Adds or updates a frequency filter
        :param qualifier:
        :param lowerbound:
        :param upperbound:
        :return:
        """
        assert qualifier is not None
        assert lowerbound is not None or upperbound is not None

        if qualifier not in self.value_filters.keys():
            self.value_filters[qualifier] = RangeCheck(lower=lowerbound, upper=upperbound)
        elif lowerbound is not None:
            self.value_filters[qualifier].lower = lowerbound
        elif upperbound is not None:
            self.value_filters[qualifier].upper = upperbound

    def add_frequency_filter(self, qualifier=None, lowerbound=None, upperbound=None):
        """
        Adds or updates a frequency filter
        :param qualifier:
        :param lowerbound:
        :param upperbound:
        :return:
        """
        assert qualifier is not None
        assert lowerbound is not None or upperbound is not None

        if qualifier not in self.frequency_filters.keys():
            self.frequency_filters[qualifier] = RangeCheck(lower=lowerbound, upper=upperbound)
        elif lowerbound is not None:
            self.frequency_filters[qualifier].lower = lowerbound
        elif upperbound is not None:
            self.frequency_filters[qualifier].upper = upperbound

    def __str__(self):
        s = "SampleFilter("
        for key, check in self.value_filters.items():
            s += "value " + str(key) + " " + str(check)
        for key, check in self.frequency_filters.items():
            s += "value " + str(key) + " " + str(check)
        s += ")"
        return s


def apply_filters(provider=None, filters=None):
    """
    Apply the specified filters and return the relevant data.   

    :param provider: the variant provider
    :param filters: filter to apply to the variant
    :return: the filtered variant
    """
    assert filters is not None
    for v in provider:
        tmp = v
        for f in filters:
            tmp = f(tmp)
            if tmp is None:
                break
        if tmp is not None and tmp.statistics.number_of_alleles > 0:
            yield tmp


if __name__ == '__main__':
    pass


