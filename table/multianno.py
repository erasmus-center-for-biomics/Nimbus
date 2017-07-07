import re
import table


class MultiAnno(table.Table):
    """A parser for MultiAnno files."""

    def __init__(self, *args, **kwargs):
        super(MultiAnno, self).__init__(*args, **kwargs)
        self.header = []

    def next(self):
        """
        Get the next entry from the table
        """
        values = super(MultiAnno, self).next()
        if len(self.header) == 0:
            try:
                # line 1
                idx = values.index("Otherinfo")
                self.header = values[0:idx]

                # line 2
                values = super(MultiAnno, self).next()
                self.header.extend(values[idx:])
            except LookupError:
                raise ValueError("Input file is not a valid multianno file")

            # get data
            values = super(MultiAnno, self).next()

        # return the data with the fields zipped
        values = [table.try_number(x) for x in values]
        return zip(self.header, values)
