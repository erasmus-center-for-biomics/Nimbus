import filters


class EqualsFilter(filters.FilterBase):

    def __init__(self, value=None):
        super(EqualsFilter, self).__init__()
        self._equals = value

    def check(self, *args, **kwargs):
        if self._equals is not None:
            if str(args[0]) != self._equals:
                raise filters.RejectedException()

class NotEqualsFilter(filters.FilterBase):

    def __init__(self, value=None):
        super(NotEqualsFilter, self).__init__()
        self._equals = value

    def check(self, *args, **kwargs):
        if self._equals is not None:
            if str(args[0]) == self._equals:
                raise filters.RejectedException()
