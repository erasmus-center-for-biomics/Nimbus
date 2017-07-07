import filters
import re


class RegExpFilter(filters.FilterBase):

    def __init__(self, value=None):
        super(RegExpFilter, self).__init__()
        self._matcher = re.compile(value)

    def check(self, *args, **kwargs):
        if not self._matcher.search(str(args[0])):
            raise filters.RejectedException()

class NotRegExpFilter(filters.FilterBase):

    def __init__(self, value=None):
        super(NotRegExpFilter, self).__init__()
        self._matcher = re.compile(value)

    def check(self, *args, **kwargs):
        if self._matcher.search(str(args[0])):
            raise filters.RejectedException()



