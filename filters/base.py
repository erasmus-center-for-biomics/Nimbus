

class RejectedException(Exception):
    pass


class FilterBase(object):

    def __init__(self, *args, **kwargs):
        pass

    def __call__(self, *args, **kwargs):
        return self.check(*args, **kwargs)

    def check(self, *args, **kwargs):
        raise NotImplementedError("Needs to be implemented")
