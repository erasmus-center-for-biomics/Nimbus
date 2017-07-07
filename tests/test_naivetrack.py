#!/bin/env python

import StringIO
import unittest
import naivetrack


class TestSimpleNaivetrackParser(unittest.TestCase):

    def setUp(self):
        # we will use a global chromosome list
        self.chromosome_list = []
        self.file = StringIO.StringIO("""#
#
chr1\t1\tT\t1\t1\t20
\tG\t1\t20\tA, forward
chr1\t20000\tT\t1\t1\t20
\tG\t1\t20\tA, forward
chr1\t20001\tT\t1\t1\t20
\tG\t1\t20\tA, forward
chr2\t20001\tG\t1\t1\t20
\tA\t1\t20\tA, forward
""")

    def runTest(self):
        parser = naivetrack.TrackParser(chromosome_list=self.chromosome_list)
        parser.register_handle(self.file)
        parser.set_maximum_chromosome_size(500000)

        entry = parser.next()
        self.assertEqual(entry.position.coordinate, 1)
        self.assertEqual(entry.position.sequence, 'T')
        self.assertEqual(entry.position.get_chromosome(), 'chr1')
        self.assertEqual(len(parser.buffer), 1)
        self.assertEqual(len(entry.alleles),1)

        entry = parser.next()
        self.assertEqual(entry.position.coordinate, 20000)
        self.assertEqual(entry.position.sequence, 'T')
        self.assertEqual(entry.position.get_chromosome(), 'chr1')
        self.assertEqual(len(parser.buffer), 2)
        self.assertEqual(len(entry.alleles),1)

        entry = parser.next()
        self.assertEqual(entry.position.coordinate, 20001)
        self.assertEqual(entry.position.sequence, 'T')
        self.assertEqual(entry.position.get_chromosome(), 'chr1')
        self.assertEqual(len(parser.buffer), 2)
        self.assertEqual(len(entry.alleles),1)

        entry = parser.next()
        self.assertEqual(entry.position.coordinate, 20001)
        self.assertEqual(entry.position.sequence, 'G')
        self.assertEqual(entry.position.get_chromosome(), 'chr2')
        self.assertEqual(len(parser.buffer), 1)
        self.assertEqual(len(entry.alleles),1)

        # end of file
        self.assertRaises(StopIteration, parser.next)

    def tearDown(self):
        pass


class TestFilterNaivetrackParser(unittest.TestCase):

    def setUp(self):
        # we will use a global chromosome list
        self.chromosome_list = []
        self.file = StringIO.StringIO("""#
#
chr1\t1\tT\t3\t3\t100
\tC\t1\t5\tA, forward
\tG\t1\t40\tA, forward
\tT\t1\t55\tA, forward
""")

    def runTest(self):
        parser = naivetrack.TrackParser(chromosome_list=self.chromosome_list)
        parser.register_handle(self.file)
        parser.set_maximum_chromosome_size(500000)

        sample_filter = naivetrack.SampleFilter()
        sample_filter.add_frequency_filter("__qual__", 0.1)
        provider = naivetrack.apply_filters(parser, filters=[sample_filter])

        entry = provider.next()
        self.assertEqual(entry.position.coordinate, 1)
        self.assertEqual(entry.position.sequence, 'T')
        self.assertEqual(entry.position.get_chromosome(), 'chr1')
        self.assertEqual(len(parser.buffer), 1)
        self.assertEqual(len(entry.alleles), 2)

        self.assertEqual(entry.alleles[0].sequence, 'G')
        self.assertEqual(entry.alleles[0].depth, 1)
        self.assertEqual(entry.alleles[0].quality, 40)

        self.assertEqual(entry.alleles[1].sequence, 'T')
        self.assertEqual(entry.alleles[1].depth, 1)
        self.assertEqual(entry.alleles[1].quality, 55)

        self.assertEqual(entry.statistics.number_of_alleles, 2)
        self.assertEqual(entry.statistics.depth, 2)
        self.assertEqual(entry.statistics.quality, 95)


class TestComplexNaivetrackParser(unittest.TestCase):

    def setUp(self):
        # we will use a global chromosome list
        self.chromosome_list = []
        self.file_a = StringIO.StringIO("""#
#
chr1\t1\tT\t1\t1\t20
\tG\t1\t20\tA, forward
chr1\t20000\tT\t1\t1\t20
\tG\t1\t20\tA, forward
chr1\t20001\tT\t1\t1\t20
\tG\t1\t20\tA, forward
chr2\t20001\tG\t1\t1\t20
\tA\t1\t20\tA, forward
""")
        self.file_b = StringIO.StringIO("""chr1\t3\tT\t1\t1\t20
\tG\t1\t20\tB, forward
chr1\t20000\tT\t1\t1\t20
\tG\t1\t20\tB, forward
chr1\t200000\tT\t1\t1\t20
\tG\t1\t20\tB, forward
chr2\t20001\tG\t1\t1\t20
\tA\t1\t20\tB, forward
chr3\t20001\tG\t1\t1\t20
\tA\t1\t20\tB, forward
""")

    def runTest(self):
        parser = naivetrack.TrackParser(chromosome_list=self.chromosome_list)
        parser.register_handle(self.file_a)
        parser.register_handle(self.file_b)
        parser.set_maximum_chromosome_size(500000)

        entry = parser.next()
        self.assertEqual(entry.position.coordinate, 1)
        self.assertEqual(entry.statistics.number_of_alleles, 1)
        self.assertEqual(entry.position.sequence, 'T')
        self.assertEqual(len(parser.buffer), 2)

        # second position
        entry = parser.next()
        self.assertEqual(entry.position.coordinate, 3)
        self.assertEqual(entry.statistics.number_of_alleles, 1)
        self.assertEqual(entry.position.sequence, 'T')

        # third position merged 2 entries
        entry = parser.next()
        self.assertEqual(entry.position.coordinate, 20000)
        self.assertEqual(entry.statistics.number_of_alleles, 2)
        self.assertEqual(entry.position.sequence, 'T')
        self.assertEqual(len(parser.buffer), 2)

        # fourth position
        entry = parser.next()
        self.assertEqual(entry.position.coordinate, 20001)
        self.assertEqual(entry.position.sequence, 'T')
        self.assertEqual(len(parser.buffer), 2)

        # fifth position
        entry = parser.next()
        self.assertEqual(entry.position.coordinate, 200000)
        self.assertEqual(entry.position.sequence, 'T')
        self.assertEqual(len(parser.buffer), 1)

        # sixt position: next chromosome test
        entry = parser.next()
        self.assertEqual(entry.position.coordinate, 20001)
        self.assertEqual(entry.position.sequence, 'G')
        self.assertEqual(entry.position.get_chromosome(), 'chr2')
        self.assertEqual(len(parser.buffer), 1)

        # sixt position: next chromosome test
        entry = parser.next()
        self.assertEqual(entry.position.coordinate, 20001)
        self.assertEqual(entry.position.sequence, 'G')
        self.assertEqual(entry.position.get_chromosome(), 'chr3')
        self.assertEqual(len(parser.buffer), 1)

        # end of file
        self.assertRaises(StopIteration, parser.next)


class TestNaivetrackSummary(unittest.TestCase):

    def setUp(self):
        # we will use a global chromosome list
        self.chromosome_list = []
        self.file = StringIO.StringIO("""#
#
chr1\t1\tT\t4\t34\t260
\tG\t1\t20\tX, forward
\tG\t3\t20\tX, reverse
\tT\t10\t20\tX, forward
\tG\t20\t200\tY, forward
chr1\t20000\tT\t1\t1\t20
\tG\t1\t20\tX, forward
""")

    def runTest(self):
        summarize = naivetrack.SummarizeEntries()
        summarize_with_strand = naivetrack.SummarizeEntries()
        summarize_with_strand.add_filter(2, "forward", "__quality_f__", "__depth_f__")

        parser = naivetrack.TrackParser(chromosome_list=self.chromosome_list)
        parser.register_handle(self.file)
        parser.set_maximum_chromosome_size(500000)

        entry = parser.next()
        self.assertEqual(entry.position.coordinate, 1)
        self.assertEqual(entry.position.sequence, 'T')
        self.assertEqual(entry.position.get_chromosome(), 'chr1')
        self.assertEqual(len(entry.alleles), 4)

        dsum = summarize(entry)
        self.assertEqual(dsum['T']['X']['__quality__'], 20)
        self.assertEqual(dsum['T']['X']['__depth__'], 10)
        self.assertEqual(dsum['G']['X']['__quality__'], 40)
        self.assertEqual(dsum['G']['X']['__depth__'], 4)
        self.assertEqual(dsum['G']['Y']['__quality__'], 200)
        self.assertEqual(dsum['G']['Y']['__depth__'], 20)

        strandsum = summarize_with_strand(entry)
        self.assertEqual(strandsum['T']['X']['__quality__'], 20)
        self.assertEqual(strandsum['T']['X']['__quality_f__'], 20)
        self.assertEqual(strandsum['T']['X']['__depth__'], 10)
        self.assertEqual(strandsum['T']['X']['__depth_f__'], 10)

        self.assertEqual(strandsum['G']['X']['__quality__'], 40)
        self.assertEqual(strandsum['G']['X']['__quality_f__'], 20)
        self.assertEqual(strandsum['G']['X']['__depth__'], 4)
        self.assertEqual(strandsum['G']['X']['__depth_f__'], 1)

        self.assertEqual(strandsum['G']['Y']['__quality__'], 200)
        self.assertEqual(strandsum['G']['Y']['__quality_f__'], 200)
        self.assertEqual(strandsum['G']['Y']['__depth__'], 20)
        self.assertEqual(strandsum['G']['Y']['__depth_f__'], 20)

        entry = parser.next()
        self.assertEqual(entry.position.coordinate, 20000)
        self.assertEqual(entry.position.sequence, 'T')
        self.assertEqual(entry.position.get_chromosome(), 'chr1')
        self.assertEqual(len(entry.alleles), 1)

        # end of file
        self.assertRaises(StopIteration, parser.next)

    def tearDown(self):
        pass

if __name__ == '__main__':
    unittest.main()
