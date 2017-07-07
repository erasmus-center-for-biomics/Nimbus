#!/bin/env python

import unittest
import StringIO
import table


class TestTable(unittest.TestCase):

    def setUp(self):
        self.data = StringIO.StringIO("""
header_a,header_b,header_c
line_1a,line_1b,line_1c
line_2a,line_2b,line_2c
""")

    def runTest(self):
        it = table.Table(h=self.data, sep=",")

        value = it.next()
        self.assertEqual(value, ["header_a", "header_b", "header_c"])

        value = it.next()
        self.assertEqual(value, ["line_1a", "line_1b", "line_1c"])

        value = it.next()
        self.assertEqual(value, ["line_2a", "line_2b", "line_2c"])

        with self.assertRaises(StopIteration):
            value = it.next()

class TestConstantWidthTable(unittest.TestCase):

    def setUp(self):
        self.data = StringIO.StringIO("""
header_a,header_b,header_c
line_1a,line_1b
line_2a,line_2b,line_2c
""")

    def runTest(self):
        it = table.ConstantWidthTable(h=self.data, sep=",")

        value = it.next()
        self.assertEqual(value, ["header_a", "header_b", "header_c"])

        with self.assertRaises(ValueError):
            value = it.next()
            print value

        value = it.next()
        self.assertEqual(value, ["line_2a", "line_2b", "line_2c"])




