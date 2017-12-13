import numpy as np
from util import MelodistTestCase


class TestRadiation(MelodistTestCase):
    def nighttime_radiation_iszero(self, s, tol=0):
        return np.allclose(s[(s.index.hour >= 22) | (s.index.hour <= 4)], 0, atol=tol)

    def test_pot_rad(self):
        station = self.station
        station.disaggregate_radiation('pot_rad')

        r = station.data_daily.glob
        rd = station.data_disagg.glob
        rdd = rd.resample('D').mean()

        assert np.all(np.isfinite(rd))
        assert np.allclose(r, rdd, atol=1e-3, equal_nan=True)
        assert self.nighttime_radiation_iszero(rd)

    def test_pot_rad_via_bc(self):
        station = self.station
        station.disaggregate_radiation('pot_rad_via_bc')

        rd = station.data_disagg.glob

        assert np.all(np.isfinite(rd))
        assert self.nighttime_radiation_iszero(rd)

    def test_mean_course(self):
        station = self.station
        station.statistics.calc_radiation_stats()
        station.disaggregate_radiation('mean_course')

        r = station.data_daily.glob
        rd = station.data_disagg.glob

        assert np.all(np.isfinite(rd))
        assert self.nighttime_radiation_iszero(rd, tol=10)
        assert np.allclose(r.resample('M').mean(), rd.resample('M').mean(), rtol=0.1)
