import tempfile

from numpy.testing import assert_equal
from pandas.testing import assert_frame_equal, assert_series_equal
from util import MelodistTestCase

import melodist


class TestIO(MelodistTestCase):
    def test_json(self):
            ss = self.station.statistics
            ss.calc_temperature_stats()
            ss.calc_precipitation_stats()
            ss.calc_humidity_stats()
            ss.calc_radiation_stats()
            ss.calc_wind_stats()

            with tempfile.NamedTemporaryFile() as tmp:
                ss.to_json(tmp.name)
                tmp.seek(0)
                ss2 = melodist.StationStatistics.from_json(tmp.name)

            assert_series_equal_kwargs = dict(
                check_dtype=False,
                check_index_type=False,
            )
            assert_frame_equal_kwargs = dict(
                check_index_type=False,
                check_column_type=False,
            )

            assert_series_equal(
                ss.temp.max_delta,
                ss2.temp.max_delta,
                **assert_series_equal_kwargs,
            )
            assert_frame_equal(
                ss.temp.mean_course,
                ss2.temp.mean_course,
                **assert_frame_equal_kwargs,
            )

            assert_equal(ss.precip.months, ss2.precip.months)
            assert all([cs1 == cs2 for cs1, cs2 in zip(ss.precip.stats, ss2.precip.stats)])

            assert ss.hum.a0 == ss2.hum.a0
            assert ss.hum.a1 == ss2.hum.a1
            assert ss.hum.kr == ss2.hum.kr
            assert_series_equal(
                ss.hum.month_hour_precip_mean,
                ss2.hum.month_hour_precip_mean,
                **assert_series_equal_kwargs,
            )

            assert_frame_equal(ss.glob.angstroem, ss2.glob.angstroem, **assert_frame_equal_kwargs)
            assert_frame_equal(ss.glob.bristcamp, ss2.glob.bristcamp, **assert_frame_equal_kwargs)
            assert_frame_equal(
                ss.glob.mean_course,
                ss2.glob.mean_course,
                **assert_frame_equal_kwargs,
            )

            assert ss.wind.a == ss2.wind.a
            assert ss.wind.b == ss2.wind.b
            assert ss.wind.t_shift == ss2.wind.t_shift
