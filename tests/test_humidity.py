import numpy as np
from util import MelodistTestCase


class TestHumidity(MelodistTestCase):
    def test_generic(self):
        for method in ('equal', 'minimal', 'dewpoint_regression',
                       'linear_dewpoint_variation', 'min_max', 'month_hour_precip_mean'):
            station = self.station
            station.statistics.calc_humidity_stats()
            station.disaggregate_temperature()
            station.disaggregate_humidity(method)

            hd = station.data_disagg.hum

            assert hd.notnull().sum() > 0
            assert not np.any(hd <= 0)
            assert not np.any(hd > 100)

    def test_equal(self):
        station = self.station
        station.disaggregate_humidity('equal')

        hd = station.data_disagg.hum
        hdd = hd.resample('D').mean()

        assert np.allclose(
            hdd,
            station.statistics.data.hum.resample('D').mean(),
            atol=1e-3,
            equal_nan=True,
        )
