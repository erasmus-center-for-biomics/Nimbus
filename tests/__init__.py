import unittest
import sys
import os.path

print(sys.version)
print(sys.executable)

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "../"))
from tests.test_naivetrack import TestSimpleNaivetrackParser, TestFilterNaivetrackParser, TestNaivetrackSummary
from tests.test_naivetrack_methods import TestNaivetrackSummaries
from tests.test_table import TestTable, TestConstantWidthTable
from tests.test_multianno import TestMultiAnno
from tests.test_filters import TestFilters
from tests.test_vcf import TestVCF
from tests.test_location import TestLocationWriter, TestLocationReader
from tests.test_variant_writer import TestVariantWriter

if __name__ == "__main__":
    unittest.main(verbosity=2)
