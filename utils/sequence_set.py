import bisect


class SequentialSet(set):
    '''
    Enabled binary search within error window
    '''

    def __init__(self, a):
        if isinstance(a, list):
            self.sset = list(set(a))
        elif isinstance(a, set):
            self.sset = list(a)
        else:
            self.sset = a

        self.sset.sort()
        self.set_relative_error()

    def set_abs_error(self, error):
        self.get_error = lambda x: error

    def set_relative_error(self, error=20e-6):
        self.get_error = lambda x: abs(x) * error

    def __contains__(self, x):
        # https://docs.python.org/3/library/bisect.html
        # Find the first value greater than or equals to the left interval. Index might right-overflow.
        left = bisect.bisect_left(self.sset, x - self.get_error(x))
        # Find the last value less than or equals to the right interval. Index might left-overflow.
        right = bisect.bisect_left(self.sset, x + self.get_error(x))-1
        return right - left >= 0

    def __iter__(self):
        return self.sset.__iter__()

    def __str__(self):
        return self.sset.__str__()

    def __repr__(self):
        return self.sset.__repr__()
