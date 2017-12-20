import numpy as np
import util


class TestTemperature(util.MelodistTestCase):
    def test_sine_min_max(self):
        station = self.station
        station.disaggregate_temperature(method='sine_min_max')

        dd = station.data_daily
        pos = util.notna_temp_days(dd, 'minmax')
        td = util.extract_days(station.data_disagg.temp, dd.index[pos])
        tdd = td.resample('D').mean().loc[pos]

        assert np.all(np.isfinite(td))
        assert np.all(tdd >= dd.tmin.min())
        assert np.all(tdd <= dd.tmax.max())
        assert (dd.temp[pos] - tdd).abs().quantile(.99) < 2

    def test_mean_course_min_max(self):
        station = self.station
        station.statistics.calc_temperature_stats()
        station.disaggregate_temperature(method='mean_course_min_max')

        dd = station.data_daily
        pos = util.notna_temp_days(dd, 'minmax')
        td = util.extract_days(station.data_disagg.temp, dd.index[pos])
        tdd = td.resample('D').mean().loc[pos]
        tdd_min = td.resample('D').min().loc[pos]
        tdd_max = td.resample('D').max().loc[pos]

        assert np.all(np.isfinite(td))
        assert np.all(tdd >= dd.tmin.min())
        assert np.all(tdd <= dd.tmax.max())
        assert np.allclose(dd.tmin[pos], tdd_min, atol=1e-3)
        assert np.allclose(dd.tmax[pos], tdd_max, atol=1e-3)

    def test_sine_mean(self):
        station = self.station
        station.disaggregate_temperature(method='sine_mean')

        dd = station.data_daily
        pos = util.notna_temp_days(dd, 'mean')
        td = util.extract_days(station.data_disagg.temp, dd.index[pos])
        tdd = td.resample('D').mean().loc[pos]
        tdd_min = td.resample('D').min().loc[pos]
        tdd_max = td.resample('D').max().loc[pos]

        assert np.all(np.isfinite(td))
        assert np.allclose(dd.temp[pos], tdd, atol=1e-3)
        assert np.allclose(
            (dd.tmax - dd.tmin)[pos],
            tdd_max - tdd_min,
            atol=1e-3
        )

    def test_mean_course_mean(self):
        station = self.station
        station.statistics.calc_temperature_stats()
        station.disaggregate_temperature(method='mean_course_mean')

        dd = station.data_daily
        pos = util.notna_temp_days(dd, 'mean')
        td = util.extract_days(station.data_disagg.temp, dd.index[pos])
        tdd = td.resample('D').mean().loc[pos]
        tdd_min = td.resample('D').min().loc[pos]
        tdd_max = td.resample('D').max().loc[pos]

        assert np.all(np.isfinite(td))
        assert np.allclose(dd.temp[pos], tdd, atol=1e-3)
        assert np.allclose(
            (dd.tmax - dd.tmin)[pos],
            tdd_max - tdd_min,
            atol=1e-3
        )
