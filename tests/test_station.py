import copy

import numpy as np
from util import MelodistTestCase, setup_station


class TestStation(MelodistTestCase):
    def setUp(self):
        self._station = setup_station()

    @property
    def station(self):
        return copy.copy(self._station)

    def test_copy(self):
        s1 = self.station
        s2 = copy.deepcopy(s1)

        s1.statistics.temp.max_delta = 3
        s1.data_disagg.temp[:] = 273.15

        assert s1.statistics.temp.max_delta != s2.statistics.temp.max_delta
        assert not np.allclose(s1.data_disagg, s2.data_disagg, equal_nan=True)
