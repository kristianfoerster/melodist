# -*- coding: utf-8 -*-
###############################################################################################################
# This file is part of MELODIST - MEteoroLOgical observation time series DISaggregation Tool                  #
# a program to disaggregate daily values of meteorological variables to hourly values                         #
#                                                                                                             #
# Copyright (C) 2016  Florian Hanzer (1,2), Kristian Förster (1,2), Benjamin Winter (1,2), Thomas Marke (1)   #
#                                                                                                             #
# (1) Institute of Geography, University of Innsbruck, Austria                                                #
# (2) alpS - Centre for Climate Change Adaptation, Innsbruck, Austria                                         #
#                                                                                                             #
# MELODIST is free software: you can redistribute it and/or modify                                            #
# it under the terms of the GNU General Public License as published by                                        #
# the Free Software Foundation, either version 3 of the License, or                                           #
# (at your option) any later version.                                                                         #
#                                                                                                             #
# MELODIST is distributed in the hope that it will be useful,                                                 #
# but WITHOUT ANY WARRANTY; without even the implied warranty of                                              #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the                                                #
# GNU General Public License for more details.                                                                #
#                                                                                                             #
# You should have received a copy of the GNU General Public License                                           #
# along with this program.  If not, see <http://www.gnu.org/licenses/>.                                       #
#                                                                                                             #
###############################################################################################################

from __future__ import print_function, division, absolute_import
import melodist
from melodist.util.bunch import Bunch
import json
import numpy as np
import pandas as pd


class StationStatistics(object):
    """
    Class representing objects that include statistical information about
    diurnal features of selected variables. A StationStatistics object is
    generally associated to a Station object for which this infomration
    is valid.
    """
    def __init__(self, data=None, lon=None, lat=None, timezone=None):
        self._data = None
        self._lon = lon
        self._lat = lat
        self._timezone = timezone

        if data is not None:
            self.data = data

        self.wind = Bunch(a=None, b=None, t_shift=None)
        self.precip = Bunch(months=None, stats=None)
        self.hum = Bunch(a0=None, a1=None, kr=None, month_hour_precip_mean=None)
        self.temp = Bunch(max_delta=None, mean_course=None)

        angstroem_df = pd.DataFrame(index=np.arange(12) + 1, columns=['a', 'b'])
        angstroem_df.a = 0.25
        angstroem_df.b = 0.75

        bristcamp_df = pd.DataFrame(index=np.arange(12) + 1, columns=['a', 'c'])
        bristcamp_df.a = 0.75
        bristcamp_df.c = 2.4

        self.glob = Bunch(angstroem=angstroem_df, bristcamp=bristcamp_df, mean_course=None)

    @property
    def data(self):
        return self._data

    @data.setter
    def data(self, df):
        assert isinstance(df, pd.DataFrame)
        assert df.index.is_all_dates
        assert df.index.resolution == 'hour'

        self._data = df.copy()

    def calc_precipitation_stats(self, months=None, avg_stats=True, percentile=50):
        """
        Calculates precipitation statistics for the cascade model while aggregating hourly observations

        Parameters
        ----------
        months :        Months for each seasons to be used for statistics (array of numpy array, default=1-12, e.g., [np.arange(12) + 1])
        avg_stats :     average statistics for all levels True/False (default=True)
        percentile :    percentil for splitting the dataset in small and high intensities (default=50)
        
        """
        if months is None:
            months = [np.arange(12) + 1]

        self.precip.months = months
        self.precip.stats = melodist.build_casc(self.data, months=months, avg_stats=avg_stats, percentile=percentile)

    def calc_wind_stats(self):
        """
        Calculates statistics in order to derive diurnal patterns of wind speed
        """
        a, b, t_shift = melodist.fit_cosine_function(self.data.wind)
        self.wind.update(a=a, b=b, t_shift=t_shift)

    def calc_humidity_stats(self):
        """
        Calculates statistics in order to derive diurnal patterns of relative humidity.
        """
        a1, a0 = melodist.calculate_dewpoint_regression(self.data, return_stats=False)
        self.hum.update(a0=a0, a1=a1)
        self.hum.kr = 12

        self.hum.month_hour_precip_mean = melodist.calculate_month_hour_precip_mean(self.data)

    def calc_temperature_stats(self):
        """
        Calculates statistics in order to derive diurnal patterns of temperature
        """
        self.temp.max_delta = melodist.get_shift_by_data(self.data.temp, self._lon, self._lat, self._timezone)
        self.temp.mean_course = melodist.util.calculate_mean_daily_course_by_month(self.data.temp, normalize=True)

    def calc_radiation_stats(self, data_daily=None, day_length=None, how='all'):
        """
        Calculates statistics in order to derive solar radiation from sunshine duration or
        minimum/maximum temperature.

        Parameters
        ----------
        data_daily : DataFrame, optional
            Daily data from the associated ``Station`` object.

        day_length : Series, optional
            Day lengths as calculated by ``calc_sun_times``.
        """
        assert how in ('all', 'seasonal', 'monthly')

        self.glob.mean_course = melodist.util.calculate_mean_daily_course_by_month(self.data.glob)

        if data_daily is not None:
            pot_rad = melodist.potential_radiation(
                melodist.util.hourly_index(data_daily.index),
                self._lon, self._lat, self._timezone)
            pot_rad_daily = pot_rad.resample('D').mean()
            obs_rad_daily = self.data.glob.resample('D').mean()

            if how == 'all':
                month_ranges = [np.arange(12) + 1]
            elif how == 'seasonal':
                month_ranges = [[3, 4, 5], [6, 7, 8], [9, 10, 11], [12, 1, 2]]
            elif how == 'monthly':
                month_ranges = zip(np.arange(12) + 1)

            def myisin(s, v):
                return pd.Series(s).isin(v).values

            def extract_months(s, months):
                return s[myisin(s.index.month, months)]

            if 'ssd' in data_daily and day_length is not None:
                for months in month_ranges:
                    a, b = melodist.fit_angstroem_params(
                        extract_months(data_daily.ssd, months),
                        extract_months(day_length, months),
                        extract_months(pot_rad_daily, months),
                        extract_months(obs_rad_daily, months),
                    )

                    for month in months:
                        self.glob.angstroem.loc[month] = a, b

            if 'tmin' in data_daily and 'tmax' in data_daily:
                df = pd.DataFrame(
                    data=dict(
                        tmin=data_daily.tmin,
                        tmax=data_daily.tmax,
                        pot_rad=pot_rad_daily,
                        obs_rad=obs_rad_daily,
                    )
                ).dropna(how='any')

                for months in month_ranges:
                    a, c = melodist.fit_bristow_campbell_params(
                        extract_months(df.tmin, months),
                        extract_months(df.tmax, months),
                        extract_months(df.pot_rad, months),
                        extract_months(df.obs_rad, months),
                    )

                    for month in months:
                        self.glob.bristcamp.loc[month] = a, c

    def to_json(self, filename=None):
        """
        Exports statistical data to a JSON formatted file

        Parameters
        ----------
        filename:    output file that holds statistics data
        """
        def json_encoder(obj):
            if isinstance(obj, pd.DataFrame) or isinstance(obj, pd.Series):
                if isinstance(obj.index, pd.MultiIndex):
                    obj = obj.reset_index()  # convert MultiIndex to columns

                return json.loads(obj.to_json(date_format='iso'))
            elif isinstance(obj, melodist.cascade.CascadeStatistics):
                return obj.__dict__
            elif isinstance(obj, np.ndarray):
                return obj.tolist()
            else:
                raise TypeError('%s not supported' % type(obj))

        d = dict(
            temp=self.temp,
            wind=self.wind,
            precip=self.precip,
            hum=self.hum,
            glob=self.glob
        )

        j = json.dumps(d, default=json_encoder, indent=4)

        if filename is None:
            return j
        else:
            with open(filename, 'w') as f:
                f.write(j)

    @classmethod
    def from_json(cls, filename):
        """
        Imports statistical data from a JSON formatted file

        Parameters
        ----------
        filename:    input file that holds statistics data
        """
        def json_decoder(d):
            if 'p01' in d and 'pxx' in d:  # we assume this is a CascadeStatistics object
                return melodist.cascade.CascadeStatistics.from_dict(d)

            return d

        with open(filename) as f:
            d = json.load(f, object_hook=json_decoder)

        stats = cls()

        stats.temp.update(d['temp'])
        stats.hum.update(d['hum'])
        stats.precip.update(d['precip'])
        stats.wind.update(d['wind'])
        stats.glob.update(d['glob'])

        if stats.temp.max_delta is not None:
            stats.temp.max_delta = pd.read_json(json.dumps(stats.temp.max_delta), typ='series').sort_index()

        if stats.temp.mean_course is not None:
            mc = pd.read_json(json.dumps(stats.temp.mean_course), typ='frame').sort_index()[np.arange(1, 12 + 1)]
            stats.temp.mean_course = mc.sort_index()[np.arange(1, 12 + 1)]

        if stats.hum.month_hour_precip_mean is not None:
            mhpm = pd.read_json(json.dumps(stats.hum.month_hour_precip_mean), typ='frame').sort_index()
            mhpm = mhpm.set_index(['level_0', 'level_1', 'level_2'])  # convert to MultiIndex
            mhpm = mhpm.squeeze()  # convert to Series
            mhpm = mhpm.rename_axis([None, None, None])  # remove index labels
            stats.hum.month_hour_precip_mean = mhpm

        for var in ('angstroem', 'bristcamp', 'mean_course'):
            if stats.glob[var] is not None:
                stats.glob[var] = pd.read_json(json.dumps(stats.glob[var])).sort_index()

        if stats.glob.mean_course is not None:
            stats.glob.mean_course = stats.glob.mean_course[np.arange(1, 12 + 1)]

        return stats
