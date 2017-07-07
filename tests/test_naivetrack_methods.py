#!/bin/env python

import StringIO
import unittest
import naivetrack
import naivetrack.utils


class TestNaivetrackSummaries(unittest.TestCase):

    def setUp(self):
        # we will use a global chromosome list
        self.chromosome_list = []
        self.file = StringIO.StringIO("""#
#
chr1\t2\tT\t2\t10\t50
\tG\t4\t20\tsample_1, forward
\tT\t6\t30\tsample_1, forward
chr1\t20000\tT\t3\t40\t130
\tA\t20\t60\tsample_1, forward
\tA\t10\t30\tsample_1, reverse
\tG\t10\t40\tsample_1, forward
""")

    def runTest(self):
        parser = naivetrack.TrackParser(chromosome_list=self.chromosome_list)
        parser.register_handle(self.file)
        parser.set_maximum_chromosome_size(500000)

        entry = parser.next()
        rval = naivetrack.utils.aggregate_allele_information(
            entry, attributes=[
                naivetrack.models.AlleleFactors["Sample"],
                naivetrack.models.AlleleFactors["Allele"]
                ])

        self.assertIn("sample_1", rval.keys())

        self.assertIn("T", rval["sample_1"].keys())
        self.assertIn("G", rval["sample_1"].keys())

        self.assertIn("__qual__", rval["sample_1"]["T"].keys())
        self.assertIn("__depth__", rval["sample_1"]["T"].keys())
        self.assertIn("__n__", rval["sample_1"]["T"].keys())

        self.assertIn("__qual__", rval["sample_1"]["G"].keys())
        self.assertIn("__depth__", rval["sample_1"]["G"].keys())
        self.assertIn("__n__", rval["sample_1"]["G"].keys())

        self.assertEqual(rval["sample_1"]["T"]["__qual__"], 30)
        self.assertEqual(rval["sample_1"]["G"]["__qual__"], 20)
        self.assertEqual(rval["sample_1"]["T"]["__depth__"], 6)
        self.assertEqual(rval["sample_1"]["G"]["__depth__"], 4)
        self.assertEqual(rval["sample_1"]["T"]["__n__"], 1)
        self.assertEqual(rval["sample_1"]["G"]["__n__"], 1)

        summarize = naivetrack.SummarizeEntries()
        rval = summarize(parser.next())
        self.assertIn("A", rval.keys())
        self.assertIn("G", rval.keys())

        self.assertEqual(rval["A"]["sample_1"]["__quality__"], 90)
        self.assertEqual(rval["A"]["sample_1"]["__depth__"], 30)
        # print rval

    def tearDown(self):
        pass

