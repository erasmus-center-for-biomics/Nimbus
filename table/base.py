import sys


def try_number(value):
    """Try converting the value to a number."""
    try:
        return int(value)
    except ValueError:
        try:
            return float(value)
        except ValueError:
            return value


class Table(object):
    """Reads an Annovar multianno output file."""

    def __init__(self, h=sys.stdin, sep="\t", *args, **kwargs):
        """Initialize the new object."""
        self._in = h
        self._sep = sep
        self._linecount = 0

    def __del__(self):
        if self._in is not sys.stdin:
            if not self._in.closed:
                self._in.close()

    def __iter__(self):
        return self

    def __next__(self):
        return self.next()

    def next(self):
        # most iterations should be a single line
        while True:
            line = self._in.next()
            line = line.rstrip('\n\r')
            self._linecount += 1

            # skip empty lines
            if len(line) == 0:
                continue

            # yield the columns
            return line.split(self._sep)


class ConstantWidthTable(Table):

    """
    A table class that checks whether a constant width is returned
    """

    def __init__(self, *args, **kwargs):
        super(ConstantWidthTable, self).__init__(*args, **kwargs)
        self._n_columns = None

    def next(self):
        values = super(ConstantWidthTable, self).next()
        if self._n_columns is None:
            self._n_columns = len(values)
        else:
            if self._n_columns != len(values):
                raise ValueError("Incorrect number of columns at line %d ([%s])" % (self._linecount, ','.join(values)))
        return values


class HeaderTable(ConstantWidthTable):

    def __init__(self, *args, **kwargs):
        super(HeaderTable, self).__init__(*args, **kwargs)
        self.header = []

    def next(self):
        values = super(HeaderTable, self).next()
        if len(self.header) == 0:
            self.header = values
            values = super(HeaderTable, self).next()
        return zip(self.header, values)


class TypedTable(HeaderTable):
    def __init__(self, *args, **kwargs):
        super(TypedTable, self).__init__(*args, **kwargs)

    def next(self):
        values = super(TypedTable, self).next()
        values = [(a, try_number(b)) for (a, b) in values]
        return values

