import filters


class SetsFilter(filters.FilterBase):

    def __init__(self, values=None):
        super(SetsFilter, self).__init__()
        self._in = values

    def check(self, *args, **kwargs):
        if self._in is not None:
            if str(args[0]) not in self._in:
                raise filters.RejectedException()
