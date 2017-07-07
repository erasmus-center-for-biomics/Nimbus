
import filters


class NumericFilter(filters.FilterBase):

    def __init__(self, minimum=None, maximum=None):
        super(NumericFilter, self).__init__()
        self._minimum = minimum
        self._maximum = maximum

    def check(self, *args, **kwargs):
        try: 
            if self._minimum is not None:
                if float(args[0]) < self._minimum:
                    raise filters.RejectedException()
            if self._maximum is not None:
                if float(args[0]) > self._maximum:
                    raise filters.RejectedException()
        except ValueError:
            raise filters.RejectedException()