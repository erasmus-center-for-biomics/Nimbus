import sys
import unittest
import StringIO
import table


class TestLocationWriter(unittest.TestCase):

    def setUp(self):
        self.chromosomes = ["chr1", "chr2"]
        self.location_a = table.Loc(0, 10001, 10001)
        self.location_b = table.Loc(0, 20001, 20011)
        self.location_c = table.Loc(1, 10001, 10001)
        self.data_a = {'a': 1.0, 'b':2.0}
        self.data_b = {'a': 2.0, 'b':3.0}
        self.data_c = {'a': 1.0, '_':4.0}

    def runTest(self):
        fout = StringIO.StringIO()

        writer = table.LocationWriter(fout=fout, chromosomes=self.chromosomes)

        writer(variant=self.location_a, columns=self.data_a)
        self.assertEqual(fout.getvalue(), "")
        self.assertEqual(len(writer.buffer), 1)

        # purposely feed the locations out of order
        writer(variant=self.location_c, columns=self.data_c)
        self.assertEqual(len(writer.buffer), 2)
        self.assertEqual(fout.getvalue(), "")

        writer(variant=self.location_b, columns=self.data_b)
        self.assertEqual(len(writer.buffer), 3)
        self.assertEqual(fout.getvalue(), "")

        # we need to write as __del__ closes the buffer
        writer._write_(writeall=True)
        self.assertNotEqual(fout.getvalue(), "")
        self.assertEqual(len(writer.buffer), 0)

        # check the number of lines
        content = fout.getvalue()
        lines = content.split("\n")
        self.assertEqual(len(lines), 5)
        self.assertEqual(lines[0][0], "#")

        # check the number of columns
        columns = lines[0].split("\t")
        self.assertEqual(len(columns), 6)

        # check the sorting
        columns = lines[3].split("\t")
        self.assertEqual(columns[0], "chr2")
        self.assertEqual(columns[1], "10001")


class TestLocationReader(unittest.TestCase):
    """TestLocationReader tests the location reader."""

    def setUp(self):
        """Setup prepares the test-case."""
        self.data = StringIO.StringIO("""
chromosome,start,end,header_a,header_b
chr1,1,10000,A line 1, B line 1
chr1,20000,30000,A line 2, B line 2
chr2,100,110,A line 3, B line 3
""")

    def runTest(self):

        # test the most basic functionality of
        # the LocationReader class
        reader = table.LocationReader(
            h=self.data,
            sep=",",
            label_chr="chromosome",
            label_start="start",
            label_end="end")

        # before the first call
        self.assertEqual(len(reader.chromosome_list), 0)

        # first call
        location, values = reader.next()
        self.assertEqual(len(reader.chromosome_list), 1)
        self.assertEqual(reader.chromosome_list[0], "chr1")
        self.assertEqual(len(values), 5)
        self.assertEqual(location.chrid, 0)
        self.assertEqual(location.start, 1)
        self.assertEqual(location.end, 10000)

        # second call
        location, values = reader.next()
        self.assertEqual(len(reader.chromosome_list), 1)
        self.assertEqual(reader.chromosome_list[0], "chr1")
        self.assertEqual(len(values), 5)
        self.assertEqual(location.chrid, 0)
        self.assertEqual(location.start, 20000)
        self.assertEqual(location.end, 30000)

        # third call
        location, values = reader.next()
        self.assertEqual(len(reader.chromosome_list), 2)
        self.assertEqual(reader.chromosome_list[1], "chr2")
        self.assertEqual(len(values), 5)
        self.assertEqual(location.chrid, 1)
        self.assertEqual(location.start, 100)
        self.assertEqual(location.end, 110)
