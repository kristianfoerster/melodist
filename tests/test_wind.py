import numpy as np
from util import MelodistTestCase


class TestWind(MelodistTestCase):
    def test_generic(self):
        for method in ('equal', 'cosine', 'random'):
            station = self.station
            station.statistics.calc_wind_stats()
            station.disaggregate_wind(method)

            wd = station.data_disagg.wind

            assert wd.notnull().sum() > 0
            assert not np.any(wd < 0)

            if method == 'equal':
                assert np.allclose(
                    station.data_daily.wind,
                    wd.resample('D').mean(),
                    atol=1e-3,
                    equal_nan=True,
                )
