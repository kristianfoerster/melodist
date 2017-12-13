import numpy as np
from util import MelodistTestCase


class TestTemperature(MelodistTestCase):
    def test_sine_min_max(self):
        station = self.station
        station.disaggregate_temperature(method='sine_min_max')

        td = station.data_disagg.temp
        tdd = td.resample('D').mean()

        assert np.all(np.isfinite(td))
        assert np.all(tdd >= station.data_daily.tmin)
        assert np.all(tdd <= station.data_daily.tmax)
        assert np.allclose(station.data_daily.temp, tdd, atol=2)

    def test_mean_course_min_max(self):
        station = self.station
        station.statistics.calc_temperature_stats()
        station.disaggregate_temperature(method='mean_course_min_max')

        td = station.data_disagg.temp
        tdd = td.resample('D').mean()

        assert np.all(np.isfinite(td))
        assert np.all(tdd >= station.data_daily.tmin)
        assert np.all(tdd <= station.data_daily.tmax)
        assert np.allclose(station.data_daily.tmin, td.resample('D').min(), atol=1e-2)
        assert np.allclose(station.data_daily.tmax, td.resample('D').max(), atol=1e-2)

    def test_sine_mean(self):
        station = self.station
        station.disaggregate_temperature(method='sine_mean')

        td = station.data_disagg.temp
        tdd = td.resample('D').mean()

        assert np.all(np.isfinite(td))
        assert np.allclose(station.data_daily.temp, tdd, atol=1e-2)
        assert np.allclose(
            station.data_daily.tmax - station.data_daily.tmin,
            td.resample('D').max() - td.resample('D').min(),
            atol=1e-2
        )

    def test_mean_course_mean(self):
        station = self.station
        station.statistics.calc_temperature_stats()
        station.disaggregate_temperature(method='mean_course_mean')

        td = station.data_disagg.temp
        tdd = td.resample('D').mean()

        assert np.all(np.isfinite(td))
        assert np.allclose(station.data_daily.temp, tdd, atol=1e-2)
        assert np.allclose(
            station.data_daily.tmax - station.data_daily.tmin,
            td.resample('D').max() - td.resample('D').min(),
            atol=1e-2
        )
