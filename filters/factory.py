import filters
import re


def parse_number(value):
    """
    Try converting the value to a number
    """
    try:
        return float(value)
    except ValueError:
        return None


class DictFilter(object):

    def __init__(self, key=None, fltr=None):
        assert key is not None
        assert filter is not None
        self._key = key
        self._filter = fltr

    def get_value(self, d):
        if self._key not in d.keys():
            raise LookupError
        value = d[self._key]
        if value is None:
            raise LookupError
        return value

    def check(self, d):
        value = self.get_value(d)
        self._filter(value)


class TupleListFilter(object):

    def __init__(self, key=None, fltr=None):
        assert key is not None
        assert filter is not None
        self._key = key
        self._filter = fltr

    def get_value(self, l):
        for key, value in l:
            if key == self._key:
                if value is None:
                    raise LookupError
                else:
                    return value
        raise LookupError

    def check(self, d):
        value = self.get_value(d)
        self._filter(value)


class FilterFactory(object):

    def __init__(self):
        self._regex_numeric = re.compile(r"(.+):[<]([0-9,\.,-]+)[,]([0-9,\.,-]+)[>]")
        self._regex_set = re.compile(r"(.+):\[(.+)\]")
        self._regex_regexp = re.compile(r"(.+)[~](.+)")
        self._regex_not_regexp = re.compile(r"(.+)[!][~](.+)")
        self._regex_equals = re.compile(r"(.+):(.+)")
        self._regex_not_equals = re.compile(r"(.+)!:(.+)")

    def get_filter(self, text):

        # number filter
        m = self._regex_numeric.match(text)
        if m is not None:
            flt = filters.NumericFilter(
                minimum=parse_number(m.group(2)),
                maximum=parse_number(m.group(3)))
            return (m.group(1), flt)

        # set filter
        m = self._regex_set.match(text)
        if m is not None:
            values = m.group(2).split(",")
            flt = filters.SetsFilter(values=values)
            return (m.group(1), flt)

        # not regexp filter
        m = self._regex_not_regexp.match(text)
        if m is not None:
            flt = filters.NotRegExpFilter(value=m.group(2))
            return (m.group(1), flt)

        # regexp filter
        m = self._regex_regexp.match(text)
        if m is not None:
            flt = filters.RegExpFilter(value=m.group(2))
            return (m.group(1), flt)

        # not equals filter
        m = self._regex_not_equals.match(text)
        if m is not None:
            flt = filters.NotEqualsFilter(value=m.group(2))
            return (m.group(1), flt)

        # equals filter
        m = self._regex_equals.match(text)
        if m is not None:
            flt = filters.EqualsFilter(value=m.group(2))
            return (m.group(1), flt)

        return None, None

    def get_tuplelist_filter(self, text):
        key, flt = self.get_filter(text)
        return TupleListFilter(key, flt)

    def get_dict_filter(self, text):
        key, flt = self.get_filter(text)
        return DictFilter(key, flt)
