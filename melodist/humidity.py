################################################################################
# This file is part of MELODIST - MEteoroLOgical observation time series       #
# DISaggregation Tool.                                                         #
#                                                                              #
# Copyright (C) 2016-2023 Florian Hanzer, Kristian Förster, Benjamin Winter,   #
# Thomas Marke                                                                 #
#                                                                              #
# MELODIST is free software: you can redistribute it and/or modify it under    #
# the terms of the GNU General Public License as published by the Free         #
# Software Foundation, either version 3 of the License, or (at your option)    #
# any later version.                                                           #
#                                                                              #
# MELODIST is distributed in the hope that it will be useful, but WITHOUT ANY  #
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS    #
# FOR A PARTICULAR PURPOSE. See the GNU General Public License for more        #
# details.                                                                     #
#                                                                              #
# You should have received a copy of the GNU General Public License along with #
# this program.  If not, see <http://www.gnu.org/licenses/>.                   #
################################################################################
import numpy as np
import pandas as pd

import melodist
import melodist.util.util as util


def disaggregate_humidity(
    data_daily,
    method='equal',
    temp=None,
    a0=None,
    a1=None,
    kr=None,
    month_hour_precip_mean=None,
    preserve_daily_mean=False,
):
    """general function for humidity disaggregation

    Args:
        daily_data: daily values
        method: keyword specifying the disaggregation method to be used
        temp: hourly temperature time series (necessary for some methods)
        kr: parameter for linear_dewpoint_variation method (6 or 12)
        month_hour_precip_mean: [month, hour, precip(y/n)] categorical mean values
        preserve_daily_mean: if True, correct the daily mean values of the disaggregated
            data with the observed daily means.

    Returns:
        Disaggregated hourly values of relative humidity.
    """
    assert method in (
        'equal',
        'minimal',
        'dewpoint_regression',
        'min_max',
        'linear_dewpoint_variation',
        'month_hour_precip_mean',
    ), 'Invalid option'

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
            tdew_delta = 0.5 * np.sin(
                (temp.index.hour + 1) * np.pi / kr - 3.0 * np.pi / 4.0
            )  # eq. (21) from Debele et al. (2007)

            tdew_nextday = tdew.shift(-24)
            tdew_nextday.iloc[-24:] = tdew.iloc[-24:]  # copy the last day

            # eq. (20) from Debele et al. (2007):
            # (corrected - the equation is wrong both in Debele et al. (2007) and Bregaglio et al.
            # (2010) - it should be (T_dp,day)_(d+1) - (T_dp,day)_d instead of the other way around)
            tdew += temp.index.hour / 24.0 * (tdew_nextday - tdew) + tdew_delta

        sat_vap_press_tdew = util.vapor_pressure(tdew, 100)
        sat_vap_press_t = util.vapor_pressure(temp, 100)
        hum_disagg = pd.Series(index=temp.index, data=100 * sat_vap_press_tdew / sat_vap_press_t)
    elif method == 'min_max':
        assert (
            'hum_min' in data_daily.columns and 'hum_max' in data_daily.columns
        ), 'Minimum and maximum humidity must be present in data frame'

        hmin = melodist.distribute_equally(data_daily.hum_min)
        hmax = melodist.distribute_equally(data_daily.hum_max)
        tmin = melodist.distribute_equally(data_daily.tmin)
        tmax = melodist.distribute_equally(data_daily.tmax)

        hum_disagg = hmax + (temp - tmin) / (tmax - tmin) * (hmin - hmax)
    elif method == 'month_hour_precip_mean':
        assert month_hour_precip_mean is not None

        precip_equal = melodist.distribute_equally(
            data_daily.precip
        )  # daily precipitation equally distributed to hourly values
        hum_disagg = pd.Series(index=precip_equal.index, dtype=float)
        locs = list(zip(hum_disagg.index.month, hum_disagg.index.hour, precip_equal > 0))
        hum_disagg[:] = month_hour_precip_mean.loc[locs].values

    if preserve_daily_mean:
        daily_mean_df = pd.DataFrame(
            data=dict(obs=data_daily.hum, disagg=hum_disagg.resample('D').mean())
        )
        bias = melodist.util.distribute_equally(daily_mean_df.disagg - daily_mean_df.obs)
        bias = bias.fillna(0)
        hum_disagg -= bias

    return hum_disagg.clip(0, 100)


def calculate_dewpoint_regression(hourly_data_obs, return_stats=False):
    temphum = hourly_data_obs[['temp', 'hum']]

    tdew = melodist.util.dewpoint_temperature(temphum.temp, temphum.hum).resample('D').mean()
    tmin = temphum.temp.groupby(temphum.index.normalize()).min()
    df = pd.DataFrame(data=dict(tmin=tmin, tdew=tdew)).dropna(how='any')

    return util.linregress(df.tmin, df.tdew, return_stats=return_stats)


def calculate_month_hour_precip_mean(hourly_data_obs):
    daily_precip_yesno = hourly_data_obs.precip.resample('D').sum() > 0
    daily_precip_yesno.index.name = None
    hum = hourly_data_obs.hum
    hum.index.name = None
    wet = daily_precip_yesno.loc[hum.index.date].values
    mhp_mean = hum.groupby([hum.index.month, hum.index.hour, wet]).mean()

    return mhp_mean
