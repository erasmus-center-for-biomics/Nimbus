#!/bin/env python

# library import
import operator

__author__ = "R.W.W.Brouwer"
__version__ = "1"
__copyright__ = "Copyright 2015, Erasmus MC"
__credits__ = []
__license__ = "To be decided"
__maintainer__ = "R.W.W.Brouwer"
__email__ = "r.w.w.brouwer@erasmusmc.nl"
__status__ = "testing"


def aggregate_allele_information(v=None, attributes=None):
    """
    Get the allele statistics from a variant per
      tags indexes tags. The statistics will be
      deposited in a dict at the end of the return value
      at the __n__, __qual__, and __depth__ keys.

    :param v: a tp.Variant object
    :param attributes: a list with attribute indexes to consider
    :return: a dict with allele statistics
    """

    def recursive_add(allele=None, cref=None, cattr=None):
        """
        Helper function to recursive add data to a dict

        :param allele: an allele object
        :param cref: the reference to the current level in the dict
        :param cattr: the attributes left to process
        :return:
        """
        if len(cattr) > 0:
            if allele.attribute(cattr[0]) not in cref.keys():
                cref[allele.attribute(cattr[0])] = {}
            recursive_add(allele, cref=cref[allele.attribute(cattr[0])], cattr=cattr[1:])
        else:
            if not cref:
                cref['__n__'] = 0
                cref['__depth__'] = 0
                cref['__qual__'] = 0
            cref['__n__'] += 1
            cref['__depth__'] += allele.depth
            cref['__qual__'] += allele.quality

    #
    assert v is not None
    attributes = attributes if attributes is not None else []
    rval = {}

    # recursively add all the data in the alleles to the return dict
    for a in v.alleles:
        # first divide on sequence
        recursive_add(a, rval, cattr=attributes)
    # return the dict with the aggregated values
    return rval


def aggregate_divide(a=None, b=None, trigger="__n__"):
    """
    Divide the recursive values in a by those in b

    :param a: The deepest leveled dict
    :param b: the shallower dict
    :param trigger: the key entry on which to trigger division
    :return: a dict with the levels of a with the value divided by b
    """

    def recursive_divide(ref=None, dx=None, dy=None):
        """
        A helper function to recursive divide both dictionaries

        :param ref: the reference to the current level in the result dict
        :param dx: the deepest dict
        :param dy: the shallowest dict
        :return:
        """

        # Recursively call this function untill we find the
        #  results. When we found the results add the quality and
        #  depth based frequencies.

        # have not yet reached the end of our deepest dict
        if trigger not in dx.keys():
            for k in dx.keys():

                # traverse deeper in dy if we can
                rdy = dy
                if k in dy.keys():
                    rdy = dy[k]

                # initialize a result entry if neccessary
                if k not in ref.keys():
                    ref[k] = {}

                # traverse deeper into the dicts (aka the unknown)
                recursive_divide(ref=ref[k], dx=dx[k], dy=rdy)
        else:
            for k in dx.keys():
                if k in dy.keys() and dy[k] != 0 and dy[k] is not None:
                    ref[k] = operator.truediv(dx[k], dy[k])
                else:
                    ref[k] = None

    rval = {}
    recursive_divide(rval, dx=a, dy=b)
    return rval


def aggregate_zip(a=None, b=None, value="__n__", result="__result__"):
    """
    Generates tuples with the value in a and b taken together

    :param a: the deepest dict
    :param b: the shallowest dict
    :param value: the value to zip
    :param result: the entry under which to state the result
    :return:
    """
    def recursive_zip(ref=None, dx=None, dy=None):
        """

        :param ref:
        :param dx:
        :param dy:
        :return:
        """
        if value not in dx.keys():
            for k in dx.keys():
                # traverse deeper in dy if we can
                rdy = dy
                if k in dy.keys():
                    rdy = dy[k]

                # initialize a result entry if neccessary
                if k not in ref.keys():
                    ref[k] = {}

                # traverse deeper into the dicts (aka the unknown)
                recursive_zip(ref=ref[k], dx=dx[k], dy=rdy)
        else:
            if value in dy.keys():
                ref[result] = (dx[value], dy[value])
            else:
                ref = None

    rval = {}
    recursive_zip(rval, dx=a, dy=b)
    return rval


if __name__ == '__main__':
    pass

