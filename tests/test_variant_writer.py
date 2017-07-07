import sys
import unittest
import StringIO
import table


class TestVariantWriter(unittest.TestCase):

    def setUp(self):
        self.chromosomes = ["chr1", "chr2"]
        self.variant_a = table.Variant(0, 10001, 10001, 'A', 'C')
        self.variant_b = table.Variant(0, 20001, 20011, '-', 'T')
        self.variant_c = table.Variant(1, 10001, 10001, 'A', '-')
        self.data_a = {'a': 1.0, 'b':2.0}
        self.data_b = {'a': 2.0, 'b':3.0}
        self.data_c = {'a': 1.0, '_':4.0}

    def runTest(self):
        fout = StringIO.StringIO()

        writer = table.VariantWriter(fout=fout, chromosomes=self.chromosomes)

        writer(variant=self.variant_a, columns=self.data_a)
        self.assertEqual(fout.getvalue(), "")
        self.assertEqual(len(writer.buffer), 1)

        # purposely feed the locations out of order
        writer(variant=self.variant_c, columns=self.data_c)
        self.assertEqual(len(writer.buffer), 2)
        self.assertEqual(fout.getvalue(), "")

        writer(variant=self.variant_b, columns=self.data_b)
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
        self.assertEqual(len(columns), 8)

        # check the sorting
        columns = lines[3].split("\t")
        self.assertEqual(columns[0], "chr2")
        self.assertEqual(columns[1], "10001")