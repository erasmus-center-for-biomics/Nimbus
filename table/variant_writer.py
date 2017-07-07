import functools
import operator
import table


@functools.total_ordering
class Variant(table.Loc):
    slots = ["reference", "alternate"]

    def __init__(self, chromosome, start, end, ref, alt):
        super(Variant, self).__init__(chromosome, start, end)
        self.reference = ref
        self.alternate = alt

    def __eq__(self, other):
        return (self.chrid, self.start, self.end, self.reference, self.alternate) == (other.chrid, other.start, other.end, other.reference, other.alternate)

    def __lt__(self, other):
        return (self.chrid, self.start, self.end, self.reference, self.alternate) < (other.chrid, other.start, other.end, other.reference, other.alternate)


class VariantWriter(table.LocationWriter):

    def __init__(self, fout, chromosomes):
        super(VariantWriter, self).__init__(fout, chromosomes)

    def _impute_header_(self):
        self.header = ["#chrom", "start", "end", "reference", "alternate"]
        self.header_base = 5
        self.header.extend(self._get_columns_())

    def _format_base_(self, var):
        return "%s\t%d\t%d\t%s\t%s" % (self.chromosomes[var.chrid],
                                       var.start, var.end,
                                       var.reference, var.alternate)
