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
import melodist.util
import pandas as pd
import numpy as np
import scipy.optimize

"""
This routine diaggregates daily values of global radiation data to hourly values
"""


def disaggregate_radiation(data_daily,
                           sun_times=None,
                           pot_rad=None,
                           method='pot_rad',
                           angstr_a=0.25,
                           angstr_b=0.5,
                           bristcamp_a=0.75,
                           bristcamp_c=2.4,
                           mean_course=None):
    """general function for radiation disaggregation

    Args:
        daily_data: daily values
        sun_times: daily dataframe including results of the util.sun_times function
        pot_rad: hourly dataframe including potential radiation
        method: keyword specifying the disaggregation method to be used
        angstr_a: parameter a of the Angstrom model (intercept)
        angstr_b: parameter b of the Angstrom model (slope)
        mean_course: monthly values of the mean hourly radiation course
        
    Returns:
        Disaggregated hourly values of shortwave radiation.
    """
    # check if disaggregation method has a valid value
    if method not in ('pot_rad', 'pot_rad_via_ssd', 'pot_rad_via_bc', 'mean_course'):
        raise ValueError('Invalid option')

    glob_disagg = pd.Series(index=melodist.util.hourly_index(data_daily.index), dtype=float)

    if method == 'mean_course':
        assert mean_course is not None

        pot_rad = pd.Series(index=glob_disagg.index, dtype=float)
        pot_rad[:] = mean_course.unstack().loc[list(zip(pot_rad.index.month, pot_rad.index.hour))].values
    else:
        assert pot_rad is not None

    pot_rad_daily = pot_rad.resample('D').mean()

    if method in ('pot_rad', 'mean_course'):
        globalrad = data_daily.glob
    elif method == 'pot_rad_via_ssd':
        # in this case use the Angstrom model
        globalrad = pd.Series(index=data_daily.index, data=0.)
        dates = sun_times.index[sun_times.daylength > 0]  # account for polar nights
        globalrad[dates] = angstroem(data_daily.ssd[dates], sun_times.daylength[dates],
                                     pot_rad_daily[dates], angstr_a, angstr_b)
    elif method == 'pot_rad_via_bc':
        # using data from Bristow-Campbell model
        globalrad = bristow_campbell(data_daily.tmin, data_daily.tmax, pot_rad_daily, bristcamp_a, bristcamp_c)

    globalrad_equal = globalrad.reindex(pot_rad.index, method='ffill')  # hourly values (replicate daily mean value for each hour)
    pot_rad_daily_equal = pot_rad_daily.reindex(pot_rad.index, method='ffill')
    glob_disagg = pot_rad / pot_rad_daily_equal * globalrad_equal
    glob_disagg[glob_disagg < 1e-2] = 0.

    return glob_disagg


def potential_radiation(dates, lon, lat, timezone, terrain_slope=0, terrain_slope_azimuth=0,
                        cloud_fraction=0, split=False):
    """
    Calculate potential shortwave radiation for a specific location and time.

    This routine calculates global radiation as described in:
    Liston, G. E. and Elder, K. (2006): A Meteorological Distribution System for
    High-Resolution Terrestrial Modeling (MicroMet), J. Hydrometeorol., 7, 217–234.

    Corrections for eccentricity are carried out following:
    Paltridge, G.W., Platt, C.M.R., 1976. Radiative processes in Meteorology and Climatology.
    Elsevier Scientific Publishing Company, Amsterdam, Oxford, New York.

    Parameters
    ----------
    dates : DatetimeIndex or array-like
        The dates for which potential radiation shall be calculated
    lon : float
        Longitude (degrees)
    lat : float
        Latitude (degrees)
    timezone : float
        Time zone
    terrain_slope : float, default 0
        Terrain slope as defined in Liston & Elder (2006) (eq. 12)
    terrain_slope_azimuth : float, default 0
        Terrain slope azimuth as defined in Liston & Elder (2006) (eq. 13)
    cloud_fraction : float, default 0
        Cloud fraction between 0 and 1
    split : boolean, default False
        If True, return a DataFrame containing direct and diffuse radiation,
        otherwise return a Series containing total radiation
    """
    solar_constant = 1367.
    days_per_year = 365.25
    tropic_of_cancer = np.deg2rad(23.43697)
    solstice = 173.0

    dates = pd.DatetimeIndex(dates)
    dates_hour = np.array(dates.hour)
    dates_minute = np.array(dates.minute)
    day_of_year = np.array(dates.dayofyear)

    # compute solar decline in rad
    solar_decline = tropic_of_cancer * np.cos(2.0 * np.pi * (day_of_year - solstice) / days_per_year)

    # compute the sun hour angle in rad
    standard_meridian = timezone * 15.
    delta_lat_time = (lon - standard_meridian) * 24. / 360.
    hour_angle = np.pi * (((dates_hour + dates_minute / 60. + delta_lat_time) / 12.) - 1.)

    # get solar zenith angle
    cos_solar_zenith = (np.sin(solar_decline) * np.sin(np.deg2rad(lat))
                        + np.cos(solar_decline) * np.cos(np.deg2rad(lat)) * np.cos(hour_angle))
    cos_solar_zenith = cos_solar_zenith.clip(min=0)
    solar_zenith_angle = np.arccos(cos_solar_zenith)

    # compute transmissivities for direct and diffus radiation using cloud fraction
    transmissivity_direct = (0.6 + 0.2 * cos_solar_zenith) * (1.0 - cloud_fraction)
    transmissivity_diffuse = (0.3 + 0.1 * cos_solar_zenith) * cloud_fraction

    # modify solar constant for eccentricity
    beta = 2. * np.pi * (day_of_year / days_per_year)
    radius_ratio = (1.00011 + 0.034221 * np.cos(beta) + 0.00128 * np.sin(beta)
                    + 0.000719 * np.cos(2. * beta) + 0.000077 * np.sin(2 * beta))
    solar_constant_times_radius_ratio = solar_constant * radius_ratio

    mu = np.arcsin(np.cos(solar_decline) * np.sin(hour_angle) / np.sin(solar_zenith_angle))
    cosi = (np.cos(terrain_slope) * cos_solar_zenith
            + np.sin(terrain_slope) * np.sin(solar_zenith_angle) * np.cos(mu - terrain_slope_azimuth))

    # get total shortwave radiation
    direct_radiation = solar_constant_times_radius_ratio * transmissivity_direct * cosi
    diffuse_radiation = solar_constant_times_radius_ratio * transmissivity_diffuse * cos_solar_zenith
    direct_radiation = direct_radiation.clip(min=0)

    df = pd.DataFrame(index=dates, data=dict(direct=direct_radiation, diffuse=diffuse_radiation))

    if split:
        return df
    else:
        return df.direct + df.diffuse


def bristow_campbell(tmin, tmax, pot_rad_daily, A, C):
    """calculates potential shortwave radiation based on minimum and maximum temperature

    This routine calculates global radiation as described in:
    Bristow, Keith L., and Gaylon S. Campbell: On the relationship between
    incoming solar radiation and daily maximum and minimum temperature.
    Agricultural and forest meteorology 31.2 (1984): 159-166.

    Args:
        daily_data: time series (daily data) including at least minimum and maximum temeprature
        pot_rad_daily: mean potential daily radiation
        A: parameter A of the Bristow-Campbell model
        C: parameter C of the Bristow-Campbell model
    Returns:
        series of potential shortwave radiation
    """

    assert tmin.index.equals(tmax.index)

    temp = pd.DataFrame(data=dict(tmin=tmin, tmax=tmax))
    temp = temp.reindex(pd.date_range(start=temp.index[0], end=temp.index[-1], freq='D'))
    temp['tmin_nextday'] = temp.tmin
    temp.tmin_nextday.iloc[:-1] = temp.tmin.iloc[1:].values

    temp = temp.loc[tmin.index]
    pot_rad_daily = pot_rad_daily.loc[tmin.index]

    dT = temp.tmax - (temp.tmin + temp.tmin_nextday) / 2

    dT_m_avg = dT.groupby(dT.index.month).mean()
    B = 0.036 * np.exp(-0.154 * dT_m_avg[temp.index.month])
    B.index = temp.index

    if isinstance(A, pd.Series):
        months = temp.index.month
        A = A.loc[months].values
        C = C.loc[months].values

    transmissivity = A * (1 - np.exp(-B * dT**C))
    R0 = transmissivity * pot_rad_daily

    return R0


def fit_bristow_campbell_params(tmin, tmax, pot_rad_daily, obs_rad_daily):
    """
    Fit the A and C parameters for the Bristow & Campbell (1984) model using observed daily
    minimum and maximum temperature and mean daily (e.g. aggregated from hourly values) solar
    radiation.

    Parameters
    ----------
    tmin : Series
        Observed daily minimum temperature.

    tmax : Series
        Observed daily maximum temperature.

    pot_rad_daily : Series
        Mean potential daily solar radiation.

    obs_rad_daily : Series
        Mean observed daily solar radiation.
    """
    def bc_absbias(ac):
        return np.abs(np.mean(bristow_campbell(df.tmin, df.tmax, df.pot, ac[0], ac[1]) - df.obs))

    df = pd.DataFrame(data=dict(tmin=tmin, tmax=tmax, pot=pot_rad_daily, obs=obs_rad_daily)).dropna(how='any')
    res = scipy.optimize.minimize(bc_absbias, [0.75, 2.4])  # i.e. we minimize the absolute bias

    return res.x


def angstroem(ssd, day_length, pot_rad_daily, a, b):
    """
    Calculate mean daily radiation from observed sunshine duration according to Angstroem (1924).

    Parameters
    ----------
    ssd : Series
        Observed daily sunshine duration.

    day_length : Series
        Day lengths as calculated by ``calc_sun_times``.

    pot_rad_daily : Series
        Mean potential daily solar radiation.

    a : float
        First parameter for the Angstroem model (originally 0.25).

    b : float
        Second parameter for the Angstroem model (originally 0.75).
    """
    if isinstance(a, pd.Series):
        months = ssd.index.month
        a = a.loc[months].values
        b = b.loc[months].values

    glob_day = (a + b * ssd / day_length) * pot_rad_daily

    return glob_day


def fit_angstroem_params(ssd, day_length, pot_rad_daily, obs_rad_daily):
    """
    Fit the a and b parameters for the Angstroem (1924) model using observed daily
    sunshine duration and mean daily (e.g. aggregated from hourly values) solar
    radiation.

    Parameters
    ----------
    ssd : Series
        Observed daily sunshine duration.

    day_length : Series
        Day lengths as calculated by ``calc_sun_times``.

    pot_rad_daily : Series
        Mean potential daily solar radiation.

    obs_rad_daily : Series
        Mean observed daily solar radiation.
    """
    df = pd.DataFrame(data=dict(ssd=ssd, day_length=day_length, pot=pot_rad_daily, obs=obs_rad_daily)).dropna(how='any')

    def angstroem_opt(x, a, b):
        return angstroem(x[0], x[1], x[2], a, b)

    x = np.array([df.ssd, df.day_length, df.pot])
    popt, pcov = scipy.optimize.curve_fit(angstroem_opt, x, df.obs, p0=[0.25, 0.75])

    return popt
