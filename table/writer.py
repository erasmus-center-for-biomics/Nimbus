
import sys


class Writer(object):

    """
    Writes data to a table 
    """
    def __init__(self, h=sys.stdout, columns=None, na_value="-", sep="\t"):
        self.columns = []
        if columns is not None:
            self.columns.extend(columns)
        self._out = h
        self.sep = sep
        self.na_value = na_value
        self.column_written = False
        self._linecount = 0

    def __del__(self):
        if self._out is not sys.stdout and self._out is not sys.stderr:
            if not self._out.closed:
                self._out.close()

    def _process_list(self, inlist=None):
        assert len(self.columns) == len(inlist)
        tmp = []
        for x in inlist:
            if x is not None:
                tmp.append(str(x))
            else:
                tmp.append(self.na_value)
        return tmp

    def _process_tuple_list(self, tuplelist=None):
        assert len(self.columns) == len(tuplelist)
        tmp = []
        for entry in tuplelist:
            if entry[1] is not None:
                tmp.append(str(entry[1]))
            else:
                tmp.append(self.na_value)
        return tmp

    def _process_dict(self, indict=None):
        assert len(self.columns) == len(indict.keys())
        tmp = []
        for x in self.columns:
            if indict[x] is not None:
                tmp.append(str(indict[x]))
            else:
                tmp.append(self.na_value)
        return tmp

    def add(self, inlist=None, tuplelist=None, indict=None):

        if not self.column_written:
            self.column_written = True
            self._out.write("%s\n" % self.sep.join(self.columns))

        tmp = []
        if inlist is not None:
            tmp = self._process_list(inlist)
        elif tuplelist is not None:
            tmp = self._process_tuple_list(tuplelist)
        elif indict is not None:
            tmp = self._process_dict(indict)

        if len(tmp) > 0:
            line = self.sep.join(tmp)
            self._out.write("%s\n" % line)
            self._linecount += 1
