import copy
import unittest

import pandas as pd

import melodist


class MelodistTestCase(unittest.TestCase):
    def setUp(self):
        self._station = setup_station()

    @property
    def station(self):
        return copy.copy(self._station)


def setup_station():
    df_hourly = pd.read_csv('examples/testdata.csv.gz', index_col=0, parse_dates=True)
    df_hourly = df_hourly.loc['2016-01-01':'2016-12-31']
    df_hourly.temp += 273.15

    df_daily = melodist.util.daily_from_hourly(df_hourly)

    station = melodist.Station(lon=8.86, lat=51.00, timezone=1, data_daily=df_daily)
    station.statistics = melodist.StationStatistics(data=df_hourly)

    return station


def extract_days(df, dates):
    """Extract full days from an hourly Series or DataFrame"""
    return pd.concat([df[d.strftime('%Y-%m-%d')] for d in dates])


def notna_temp_days(dd, kind):
    if kind == 'minmax':
        pos = (pd.DataFrame(data=[dd.tmin, dd.tmax,
                                  dd.tmin.shift(-1), dd.tmin.shift(1),
                                  dd.tmax.shift(-1), dd.tmax.shift(1)])
               .T
               .notna()
               .all(axis=1))
    elif kind == 'mean':
        pos = (pd.DataFrame(data=[dd.tmin, dd.tmax,
                                  dd.tmin.shift(-1), dd.tmin.shift(1),
                                  dd.tmax.shift(-1), dd.tmax.shift(1),
                                  dd.temp, dd.temp.shift(-1), dd.temp.shift(1)])
               .T
               .notna()
               .all(axis=1))

    return pos
