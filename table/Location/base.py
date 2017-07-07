#!/bin/env python

import functools
import table


@functools.total_ordering
class Loc(object):
    """An object to represent a location."""

    slots = ["chrid", "start", "end"]

    def __init__(self, chromosome, start, end):
        """New location object"""
        self.chrid = chromosome
        self.start = start
        self.end = end

    def __eq__(self, other):
        return (self.chrid, self.start, self.end) == (other.chrid, other.start, other.end)

    def __lt__(self, other):
        return (self.chrid, self.start, self.end) < (other.chrid, other.start, other.end)

    def __str__(self):
        return "chrid=%d;start=%d;end=%d" % (self.chrid, self.start, self.end)
