#
# The VCFWriter class
#

import table.VCF
import table


class VCFWriter(table.LocationWriter):
    """VCFWriter writes VCF entries to an output file."""

    def __init__(self, *args, **kwargs):
        """__init__ creates a new VCFWriter object."""
        super(VCFWriter, self).__init__(*args, **kwargs)
        self.samples = []
        self.meta = table.VCF.Meta()
        if "meta" in kwargs.keys():
            self.meta = kwargs["meta"]
        # set the header
        self.header = [
            "#CHROM",
            "POS",
            "ID",
            "REF",
            "ALT",
            "QUAL",
            "FILTER",
            "INFO"
        ]
        self.header_base = 5

    def _format_base_(self, var):
        """Format the base variant."""
        return "%s\t%d\t%s\t%s\t%s" % (
            self.chromosomes[var.chrid],
            var.start + 1,
            var.name,
            var.reference,
            ','.join(var.alternate))

    def _write_header_(self):
        """Write the header."""
        self.meta.write(self.fout)
        self.fout.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO")
        if len(self.samples) > 0:
            self.fout.write("\tFORMAT\t%s\n" % ("\t".join(self.samples)))
        self.header.append("FORMAT")
        self.header.extend(self.samples)
        self.headerwritten = True

    def __call__(self, *args, **kwargs):
        """__call__ prepares a new VCF entry for writing."""
        variant = args[0]
        qual = args[1]
        fltr = args[2]
        info = args[3]
        samples = args[4]

        # add the standard columns
        columns = []
        columns.append(("QUAL", qual))
        columns.append(("FILTER", fltr))
        columns.append(("INFO", ";".join([k + "=" + str(v) for k, v in info])))

        # add the format if there are samples
        if len(samples) > 0:
            fmt = ":".join([x for x, _ in samples[0][1]])
            columns.append(("FORMAT", fmt))

        # add the samples
        sampledict = dict(samples)
        for key in self.samples:
            columns.append((
                key,
                ":".join([v for _, v in sampledict[key]])))

        # append the data to the buffer
        self.buffer.append((variant, dict(columns)))
        if len(self.buffer) >= self.buffertrigger:
            self._write_()



