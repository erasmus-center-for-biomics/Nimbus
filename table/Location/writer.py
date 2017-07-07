import operator
import table.Location
import table


class LocationWriter(object):

    """
    Writes data to an output location
    """

    def __init__(self, fout, chromosomes, *args, **kwargs):
        self.fout = fout
        self.chromosomes = chromosomes
        self.buffer = []
        self.buffertrigger = 1000
        self.emit = 800
        self.header = []
        self.header_base = 0
        self.na = "NA"
        self.headerwritten = False
        self.impute_header = True

    def __del__(self):
        self._write_(writeall=True)

    def __call__(self, *args, **kwargs):
        """
        Adds elements to the buffer
        """
        variant = None
        columns = None
        if len(args) >= 2:
            variant = args[0]
            columns = args[1]
        else:
            variant = kwargs["variant"]
            columns = kwargs["columns"]

        self.buffer.append((variant, columns))
        if len(self.buffer) >= self.buffertrigger:
            self._write_()

    def _get_columns_(self):
        allentries = []
        for loc, cols in self.buffer:
            if cols is not None:
                allentries.extend(cols.keys())
        return sorted(list(set(allentries)))

    def _impute_header_(self):
        """
        Imputes the header from the columns in the buffer
        """
        self.header = ["#chrom", "start", "end"]
        self.header_base = 3
        self.header.extend(self._get_columns_())

    def _write_header_(self):
        if self.impute_header or len(self.header) == 0:
            self._impute_header_()

        self.fout.write("%s\n" % '\t'.join(self.header))
        self.headerwritten = True

    def _format_base_(self, var):
        return "%s\t%d\t%d" % (self.chromosomes[var.chrid], var.start, var.end)

    def _write_(self, writeall=False):
        """
        Writes data from the buffer to the output file
        """
        # write the header if not done so already
        if not self.headerwritten:
            self._write_header_()

        # sort the buffer befor each write actions
        self.buffer.sort(key=operator.itemgetter(0))

        # determine the number of entries to write
        end = self.emit if not writeall else len(self.buffer)

        # write the entries
        for i in xrange(end):
            (var, columns) = self.buffer[i]
            line = self._format_base_(var)

            # append the required columns
            extend = ""
            for k in self.header[self.header_base:]:
                if k in columns and columns[k] is not None:
                    extend += "\t%s" % str(columns[k])
                else:
                    extend += "\t%s" % self.na
            if len(extend) > 0:
                line += extend
            line += "\n"
            self.fout.write(line)

        # remove written entries from the buffer
        self.buffer = self.buffer[end:]

