
import table.VCF
import table


class VCFReader(table.Table):
    """A class to read VCF a file."""

    def __init__(self, *args, **kwargs):
        """Construct a new VCF reader object."""
        kwargs["sep"] = "\t"
        super(VCFReader, self).__init__(*args, **kwargs)
        self.meta = table.VCF.Meta()
        self.header = []
        self.samples = []
        self.chromosomes = []
        if "chromosomes" in kwargs.keys():
            self.chromosomes = kwargs["chromosomes"]

    def parse_vcfvariant(self, values):
        """Get the VCF variant."""
        # add the chromosome to the chromosome list
        try:
            chridx = self.chromosomes.index(values[0])
        except ValueError:
            # add if not found
            self.chromosomes.append(values[0])
            chridx = len(self.chromosomes) - 1

        # create the vcf variant
        var = table.VCF.Variant(
            chridx,
            int(values[1]),
            values[2],
            values[3],
            values[4].split(","))
        return var

    def parse_info(self, values):
        """Parse the info field."""
        info = []
        content = values[7].split(";")
        for elem in content:
            parts = elem.split("=")
            if len(parts) == 2:
                info.append((parts[0], parts[1]))
            else:
                info.append((parts[0], True))
        return info

    def parse_format(self, values):
        """Parse the format column if present."""
        if len(values) <= 8:
            return []
        return values[8].split(':')

    def parse_samples(self, values):
        """Parse the sample information from the VCF file."""
        retval = []
        # return an empty list if no data was present
        if len(values) <= 8:
            return retval

        # split the data on : and return a
        # double pair list where the first level
        # is the sample and the second level is
        # the format
        fvar = self.parse_format(values)
        for i in xrange(len(self.samples)):
            cidx = len(self.header) + i
            sampleid = self.samples[i]
            data = values[cidx].split(":")
            svar = zip(fvar, data)
            retval.append((sampleid, svar))
        return retval

    def interpret(self, values):
        """Interpret the values and process the results in meaningfull units."""
        vcfvar = self.parse_vcfvariant(values)
        qual = None
        try:
            qual = float(values[5])
        except ValueError:
            pass
        flt = values[6]
        info = self.parse_info(values)
        sampleinfo = self.parse_samples(values)

        # return the interpreted values
        return vcfvar, qual, flt, info, sampleinfo

    def next(self):
        """Get the next entry from the table."""
        # get the values
        values = super(VCFReader, self).next()

        # get the header and meta information if necessary
        if len(self.header) == 0:

            # parse the meta information
            while True:
                if values[0].startswith("##"):
                    self.meta.add_from_line(values[0])
                else:
                    break
                values = super(VCFReader, self).next()
            values[0] = values[0][1:]

            # get the header and refresh the values list
            self.header = values

            # if we have sample specific columns
            if len(self.header) > 8:
                self.samples = self.header[9:]
                self.header = self.header[:9]

            # set the values
            values = super(VCFReader, self).next()

        # return the interpreted data
        return self.interpret(values)
