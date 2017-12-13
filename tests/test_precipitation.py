import numpy as np
from util import MelodistTestCase


class TestPrecipitation(MelodistTestCase):
    def test_equal(self):
        station = self.station
        station.disaggregate_precipitation('equal')

        p = station.data_daily.precip
        pd = station.data_disagg.precip
        pdd = pd.resample('D').sum()

        assert np.allclose(p, pdd, atol=1e-3, equal_nan=True)

    def test_cascade(self):
        station = self.station
        station.statistics.calc_precipitation_stats()
        station.disaggregate_precipitation('cascade')

        p = station.data_daily.precip
        pd = station.data_disagg.precip
        pdd = pd.resample('D').sum()
        pos = p >= 0

        assert np.allclose(station.data_daily.precip[pos], pdd[pos], atol=1e-3, equal_nan=True)
