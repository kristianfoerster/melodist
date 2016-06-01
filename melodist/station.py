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
import melodist.util
import pandas as pd

class Station(object):
    """
    Class representing meteorological stations including all relevant 
    information such as metadata and meteorological time series (observed
    and disaggregated)
    """
    _columns_daily=[
        'tmean',
        'tmin',
        'tmax',
        'precip',
        'glob',
        'ssd',
        'hum',
        'wind',
    ]

    _columns_hourly=[
        'temp',
        'precip',
        'glob',
        'hum',
        'wind',
    ]

    def __init__(self, id=None, name=None, lon=None, lat=None, timezone=None, data_daily=None):
        self._lon = None
        self._lat = None
        self._timezone = None
        self._statistics = None
        self._data_daily = None
        self._data_disagg = None

        self.statistics = melodist.StationStatistics(lon=lon, lat=lat)
        self.id = id
        self.name = name
        self.lon = lon
        self.lat = lat
        self.timezone = timezone
        self.sun_times = None

        if data_daily is not None:
            self.data_daily = data_daily

    @property
    def data_daily(self):
        """
        Daily meteorological time series either derived through observations
        or aggregation of hourly data for testing purposes.
        """
        return self._data_daily

    @data_daily.setter
    def data_daily(self, df):
        assert isinstance(df, pd.DataFrame)
        assert df.index.is_all_dates
        # for col in df:
        #     assert col in Station._columns_daily
        assert df.index.resolution == 'day'
        assert df.index.is_monotonic_increasing

        if df.index.freq is None: # likely some days are missing
            df = df.reindex(pd.DatetimeIndex(start=df.index[0], end=df.index[-1], freq='D'))

        for var in 'tmin', 'tmax', 'tmean':
            if var in df: assert not any(df[var] < 200), 'Implausible temperature values detected - temperatures must be in K'

        self._data_daily = df.copy()

        # create data frame for disaggregated data:
        index = melodist.util.hourly_index(df.index)
        df = pd.DataFrame(index=index, columns=Station._columns_hourly, dtype=float)
        self._data_disagg = df

        if self.timezone is not None:
            self.calc_sun_times()

    @property
    def lon(self):
        """
        Longitute of the station
        """
        return self._lon

    @lon.setter
    def lon(self, lon):
        self._lon = lon
        self.statistics._lon = lon

    @property
    def lat(self):
        """
        Latitute of the station
        """
        return self._lat

    @lat.setter
    def lat(self, lat):
        self._lat = lat
        self.statistics._lat = lat

    @property
    def timezone(self):
        """
        Timezone indicates the differnce in hours calculated from UTC
        
        Negative values indicate timezones later than UTC, i.e. west of 0 deg
        long. Positive values indicate the reverse.
        """
        return self._timezone

    @timezone.setter
    def timezone(self, timezone):
        self._timezone = timezone
        self.statistics._timezone = timezone

    @property
    def statistics(self):
        """
        The associated StationStatistics object
        """

        return self._statistics

    @statistics.setter
    def statistics(self, s):
        assert isinstance(s, melodist.StationStatistics)
        s._lon = self.lon
        s._lat = self.lat
        s._timezone = self.timezone
        self._statistics = s

    @property
    def data_disagg(self):
        """
        All results derived through disaggregation will be stored in this
        property.
        """

        return self._data_disagg

    def calc_sun_times(self):
        """
        Computes the times of sunrise, solar noon, and sunset for each day. 
        """

        self.sun_times = melodist.util.get_sun_times(self.data_daily.index, self.lon, self.lat, self.timezone)


    def disaggregate_wind(self, method='equal'):
        """
        Disaggregate wind speed.

        Parameters
        ----------
        method : str, optional
            Disaggregation method.

            ``equal``
                Mean daily wind speed is duplicated for the 24 hours of the day. (Default)

            ``cosine``
                Distributes daily mean wind speed using a cosine function derived from hourly
                observations.

            ``random``
                Draws random numbers to distribute wind speed (usually not conserving the
                daily average).
        """
        self.data_disagg.wind = melodist.disaggregate_wind(self.data_daily.wind, method=method, **self.statistics.wind)

    def disaggregate_humidity(self, method='equal'):
        """
        Disaggregate relative humidity.

        Parameters
        ----------
        method : str, optional
            Disaggregation method.

            ``equal``
                Mean daily humidity is duplicated for the 24 hours of the day. (Default)

            ``minimal``:
                Calculates humidity from daily dew point temperature by setting the dew point temperature
                equal to the daily minimum temperature.

            ``dewpoint_regression``:
                Calculates humidity from daily dew point temperature by calculating dew point temperature
                using ``Tdew = a * Tmin + b``, where ``a`` and ``b`` are determined by calibration.

            ``linear_dewpoint_variation``:
                Calculates humidity from hourly dew point temperature by assuming a linear dew point
                temperature variation between consecutive days.

            ``min_max``:
                Calculates hourly humidity from observations of daily minimum and maximum humidity.
        """
        self.data_disagg.hum = melodist.disaggregate_humidity(self.data_daily, temp=self.data_disagg.temp, method=method, **self.statistics.hum)

    def disaggregate_temperature(self, method='sine', min_max_time='fix', mod_nighttime=False):
        """
        Disaggregate air temperature.

        Parameters
        ----------
        method : str, optional
            Disaggregation method.

            ``sine``
                Hourly temperatures follow a sine function. (Default)

        min_max_time : str, optional
            Method to determine the time of minimum and maximum temperature.

            ``fix``:
                Minimum/maximum temperature are assumed to occur at 07:00/14:00 local time.

            ``sun_loc``:
                Minimum/maximum temperature are assumed to occur at sunrise / solar noon + 2 h.

            ``sun_loc_shift``:
                Minimum/maximum temperature are assumed to occur at sunrise / solar noon + monthly mean shift.

        mod_nighttime : bool, optional
            Use linear interpolation between minimum and maximum temperature.
        """
        self.data_disagg.temp = melodist.disaggregate_temperature(self.data_daily, method=method, min_max_time=min_max_time, max_delta=self.statistics.temp.max_delta, sun_times=self.sun_times, mod_nighttime=mod_nighttime)

    def disaggregate_precipitation(self, method='equal', zerodiv='uniform', shift=0, master_precip=None):
        """
        Disaggregate precipitation.

        Parameters
        ----------
        method : str, optional
            Disaggregation method.

            ``equal``
                Daily precipitation is distributed equally over the 24 hours of the day. (Default)

            ``cascade``
                Hourly precipitation values are obtained using a cascade model set up using
                hourly observations.

        zerodiv : str, optional
            Method to deal with zero division, relevant for ``method='masterstation'``.

            ``uniform``
                Use uniform distribution. (Default)

        master_precip : Series, optional
            Hourly precipitation records from a representative station
            (required for ``method='masterstation'``).
        """
        if method == 'equal':
            precip_disagg = melodist.disagg_prec(self.data_daily, method=method, shift=shift)
        elif method == 'cascade':
            precip_disagg = pd.Series(index=self.data_disagg.index)

            for months, stats in zip(self.statistics.precip.months, self.statistics.precip.stats):
                precip_daily = melodist.seasonal_subset(self.data_daily.precip, months=months)
                if len(precip_daily) > 1:
                    data = melodist.disagg_prec(precip_daily, method=method, cascade_options=stats, shift=shift, zerodiv=zerodiv)
                    precip_disagg.loc[data.index] = data
        elif method == 'masterstation':
            precip_disagg = melodist.precip_master_station(self.data_daily.precip, master_precip, zerodiv)

        self.data_disagg.precip = precip_disagg

    def disaggregate_radiation(self, method='pot_rad'):
        """
        Disaggregate solar radiation.

        Parameters
        ----------
        method : str, optional
            Disaggregation method.

            ``pot_rad``
                Calculates potential clear-sky hourly radiation and scales it according to the
                mean daily radiation. (Default)

            ``pot_rad_via_ssd``
                Calculates potential clear-sky hourly radiation and scales it according to the
                observed daily sunshine duration.

            ``pot_rad_via_bc``
                Calculates potential clear-sky hourly radiation and scales it according to daily
                minimum and maximum temperature.
        """
        if self.sun_times is None:
            self.calc_sun_times()

        self.data_disagg.glob = melodist.disaggregate_radiation(
            self.data_daily,
            self.sun_times,
            melodist.potential_radiation(self.data_disagg.index, self.lon, self.lat, self.timezone),
            method=method,
            angstr_a=self.statistics.glob.angstroem_a,
            angstr_b=self.statistics.glob.angstroem_b,
            bristcamp_a=self.statistics.glob.bristcamp_a,
            bristcamp_c=self.statistics.glob.bristcamp_c
        )

    def interpolate(self, column_hours, method='linear', limit=24, limit_direction='both', **kwargs):
        """
        Wrapper function for ``pandas.Series.interpolate`` that can be used to
        "disaggregate" values using various interpolation methods.

        Parameters
        ----------
        column_hours : dict
            Dictionary containing column names in ``data_daily`` and the hour
            values they should be associated to.

        method, limit, limit_direction, **kwargs
            These parameters are passed on to ``pandas.Series.interpolate``.

        Examples
        --------
        Assume that ``mystation.data_daily.T7``, ``mystation.data_daily.T14``,
        and ``mystation.data_daily.T19`` contain air temperature measurements
        taken at 07:00, 14:00, and 19:00.
        We can use the interpolation functions provided by pandas/scipy to derive
        hourly values:

        >>> mystation.data_hourly.temp = mystation.interpolate({'T7': 7, 'T14': 14, 'T19': 19}) # linear interpolation (default)
        >>> mystation.data_hourly.temp = mystation.interpolate({'T7': 7, 'T14': 14, 'T19': 19}, method='cubic') # cubic spline
        """
        kwargs = dict(kwargs, method=method, limit=limit, limit_direction=limit_direction)
        data = melodist.util.prepare_interpolation_data(self.data_daily, column_hours)
        return data.interpolate(**kwargs)
