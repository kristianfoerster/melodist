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
import melodist.util.util as util
import numpy as np
import pandas as pd

def disaggregate_humidity(data_daily, method='equal', temp=None, a0=None, a1=None, kr=None):
    """general function for humidity disaggregation

    Args:
        daily_data: daily values
        method: keyword specifying the disaggregation method to be used
        temp: hourly temperature time series (necessary for some methods)
        kr: parameter for linear_dewpoint_variation method (6 or 12)

    Returns:
        Disaggregated hourly values of relative humidity.
    """
    assert method in ('equal', 'minimal', 'dewpoint_regression', 'min_max', 'linear_dewpoint_variation'), 'Invalid option'

    if method == 'equal':
        hum_disagg = melodist.distribute_equally(data_daily.hum)
    elif method in ('minimal', 'dewpoint_regression', 'linear_dewpoint_variation'):
        if method == 'minimal':
            a0 = 0
            a1 = 1

        assert a0 is not None and a1 is not None, 'a0 and a1 must be specified'
        tdew_daily = a0 + a1 * data_daily.tmin

        tdew = melodist.distribute_equally(tdew_daily)

        if method == 'linear_dewpoint_variation':
            assert kr is not None, 'kr must be specified'
            assert kr in (6, 12), 'kr must be 6 or 12'
            tdew_delta = 0.5 * np.sin((temp.index.hour + 1) * np.pi / kr - 3. * np.pi / 4.) # eq. (21) from Debele et al. (2007)

            tdew_nextday = tdew.shift(-24)
            tdew_nextday.iloc[-24:] = tdew.iloc[-24:] # copy the last day

            # eq. (20) from Debele et al. (2007):
            # (corrected - the equation is wrong both in Debele et al. (2007) and Bregaglio et al. (2010) - it should
            # be (T_dp,day)_(d+1) - (T_dp,day)_d instead of the other way around)
            tdew += temp.index.hour / 24. * (tdew_nextday - tdew) + tdew_delta

        sat_vap_press_tdew = util.vapor_pressure(tdew, 100)
        sat_vap_press_t = util.vapor_pressure(temp, 100)
        hum_disagg = 100 * sat_vap_press_tdew / sat_vap_press_t
    elif method == 'min_max':
        assert 'hum_min' in data_daily.columns and 'hum_max' in data_daily.columns, 'Minimum and maximum humidity must be present in data frame'

        hmin = melodist.distribute_equally(data_daily.hum_min)
        hmax = melodist.distribute_equally(data_daily.hum_max)
        tmin = melodist.distribute_equally(data_daily.tmin)
        tmax = melodist.distribute_equally(data_daily.tmax)

        hum_disagg = hmax + (temp - tmin) / (tmax - tmin) * (hmin - hmax)

    return hum_disagg.clip(0, 100)


def calculate_dewpoint_regression(hourly_data_obs, return_stats=False):
    temphum = hourly_data_obs[['temp', 'hum']]

    tdew = melodist.util.dewpoint_temperature(temphum.temp, temphum.hum).resample('D').mean()
    tmin = temphum.temp.groupby(temphum.index.date).min()
    df = pd.DataFrame(data=dict(tmin=tmin, tdew=tdew)).dropna(how='any')

    return util.linregress(df.tmin, df.tdew, return_stats=return_stats)
