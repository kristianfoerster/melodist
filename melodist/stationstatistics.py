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
        self.hum = Bunch(a0=None, a1=None, kr=None)
        self.temp = Bunch(max_delta=None)
        self.glob = Bunch(angstroem_a=0.25, angstroem_b=0.75, bristcamp_a=0.75, bristcamp_c=2.4)

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
        Calculates statistics in order to derive diurnal patterns of relative humidity (requires temperature data as well)        
        """
        a1, a0 = melodist.calculate_dewpoint_regression(self.data, return_stats=False)
        self.hum.update(a0=a0, a1=a1)
        self.hum.kr = 12

    def calc_temperature_stats(self):
        """
        Calculates statistics in order to derive diurnal patterns of temperature
        """
        self.temp.max_delta = melodist.get_shift_by_data(self.data.temp, self._lon, self._lat, self._timezone)

    def calc_radiation_stats(self, data_daily, day_length=None):
        """
        Calculates statistics in order to derive solar radiation from sunshine duration or
        minimum/maximum temperature.

        Parameters
        ----------
        data_daily : DataFrame
            Daily data from the associated ``Station`` object.

        day_length : Series, optional
            Day lengths as calculated by ``calc_sun_times``.
        """
        if 'ssd' in data_daily and day_length is not None:
            df = pd.DataFrame(data=dict(ssd=data_daily.ssd, day_length=day_length)).dropna(how='any')
            pot_rad = melodist.potential_radiation(melodist.util.hourly_index(df.index), self._lon, self._lat, self._timezone)
            pot_rad_daily = pot_rad.resample('D').mean()
            obs_rad_daily = self.data.glob.resample('D').mean()
            a, b = melodist.fit_angstroem_params(data_daily.ssd, day_length, pot_rad_daily, obs_rad_daily)
            self.glob.angstroem_a = a
            self.glob.angstroem_b = b

        if 'tmin' in data_daily and 'tmax' in data_daily:
            pot_rad = melodist.potential_radiation(melodist.util.hourly_index(df.index), self._lon, self._lat, self._timezone)
            pot_rad_daily = pot_rad.resample('D').mean()
            obs_rad_daily = self.data.glob.resample('D').mean()
            df = pd.DataFrame(
                data=dict(
                    tmin=data_daily.tmin,
                    tmax=data_daily.tmax,
                    pot_rad=pot_rad_daily,
                    obs_rad=obs_rad_daily,
                )
            ).dropna(how='any')
            a, c = melodist.fit_bristow_campbell_params(df.tmin, df.tmax, df.pot_rad, df.obs_rad)
            self.glob.bristcamp_a = a
            self.glob.bristcamp_c = c

    def to_json(self, filename=None):
        """
        Exports statistical data to a JSON formatted file

        Parameters
        ----------
        filename:    output file that holds statistics data
        """
        def json_encoder(obj):
            if isinstance(obj, pd.DataFrame) or isinstance(obj, pd.Series):
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
            if 'p01' in d and 'pxx' in d: # we assume this is a CascadeStatistics object
                return melodist.cascade.CascadeStatistics.from_dict(d)

            return d

        with open(filename) as f:
            d = json.load(f, object_hook=json_decoder)

        stats = cls()

        stats.temp.update(d['temp'])
        stats.hum.update(d['hum'])
        stats.precip.update(d['precip'])
        stats.wind.update(d['wind'])

        if stats.temp.max_delta is not None:
            stats.temp.max_delta = pd.read_json(json.dumps(stats.temp.max_delta), typ='series')

        return stats
