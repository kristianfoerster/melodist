import copy


# based on http://code.activestate.com/recipes/52308-the-simple-but-handy-collector-of-a-bunch-of-named
class Bunch(dict):
    def __init__(self, **kw):
        dict.__init__(self, kw)
        self.__dict__ = self

    def __copy__(self):
        copy = type(self)(**self)
        return copy

    def __deepcopy__(self, memo):
        return Bunch(**copy.deepcopy(dict(self)))

    def copy(self):
        return self.__copy__()
