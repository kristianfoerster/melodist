import numpy as np
import util


class TestRadiation(util.MelodistTestCase):
    def nighttime_radiation_iszero(self, s, tol=0):
        return np.allclose(s[(s.index.hour >= 22) | (s.index.hour <= 4)], 0, atol=tol)

    def test_pot_rad(self):
        station = self.station
        station.disaggregate_radiation('pot_rad')

        dd = station.data_daily
        pos = dd.glob.notna()
        r = dd.glob[pos]
        rd = util.extract_days(station.data_disagg.glob, dd.index[pos])
        rdd = rd.resample('D').mean().loc[pos]

        assert np.all(np.isfinite(rd))
        assert np.allclose(r, rdd, atol=1e-3, equal_nan=True)
        assert self.nighttime_radiation_iszero(rd)

    def test_pot_rad_via_bc(self):
        station = self.station
        station.disaggregate_radiation('pot_rad_via_bc')

        dd = station.data_daily
        pos = util.notna_temp_days(dd, 'minmax')
        rd = util.extract_days(station.data_disagg.glob, dd.index[pos])

        assert np.all(np.isfinite(rd))
        assert self.nighttime_radiation_iszero(rd)

    def test_mean_course(self):
        station = self.station
        station.statistics.calc_radiation_stats()
        station.disaggregate_radiation('mean_course')

        dd = station.data_daily
        pos = dd.glob.notna()
        r = dd.glob[pos]
        rd = util.extract_days(station.data_disagg.glob, dd.index[pos])

        assert np.all(np.isfinite(rd))
        assert self.nighttime_radiation_iszero(rd, tol=10)
        assert np.allclose(r.resample('M').mean(), rd.resample('M').mean(), rtol=0.1)
