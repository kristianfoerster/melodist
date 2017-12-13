import copy
import funcs
import unittest


class MelodistTestCase(unittest.TestCase):
    def setUp(self):
        self._station = funcs.setup_station()

    @property
    def station(self):
        return copy.copy(self._station)
