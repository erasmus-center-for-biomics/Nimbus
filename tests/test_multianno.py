import os
import os.path
import unittest
import table


class TestMultiAnno(unittest.TestCase):

    def setUp(self):
        self.infile = os.path.join(
            os.path.dirname(__file__), "../data/test.hg19_multianno.txt")
        self.handle = open(self.infile, 'r')

    def runTest(self):
        parser = table.MultiAnno(h=self.handle)

        data = parser.next()
        self.assertEqual(len(data), 96)
        self.assertEqual(data[0][0], "Chr")
        self.assertEqual(data[1][0], "Start")
        self.assertEqual(data[2][0], "End")

        self.assertEqual(data[0][1], "chr1")
        self.assertEqual(data[1][1], 13327)
        self.assertEqual(data[2][1], 13327)

        data = parser.next()
        self.assertEqual(len(data), 96)
        self.assertEqual(data[0][0], "Chr")
        self.assertEqual(data[1][0], "Start")
        self.assertEqual(data[2][0], "End")

        data = parser.next()
        self.assertEqual(len(data), 96)
        self.assertEqual(data[0][0], "Chr")
        self.assertEqual(data[1][0], "Start")
        self.assertEqual(data[2][0], "End")


    def tearDown(self):
        pass

