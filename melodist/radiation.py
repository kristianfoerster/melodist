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
from math import cos, sin, acos, pi
import melodist.util
import pandas as pd
import numpy as np
import scipy.optimize

"""
This routine diaggregates daily values of global radiation data to hourly values
"""


def disaggregate_radiation(data_daily, sun_times, pot_rad, method='pot_rad', angstr_a=0.25, angstr_b=0.5, bristcamp_a=0.75, bristcamp_c=2.4):
    """general function for radiation disaggregation

    Args:
        daily_data: daily values
        sun_times: daily dataframe including results of the util.sun_times function
        pot_rad: hourly dataframe including potential radiation
        method: keyword specifying the disaggregation method to be used
        angstr_a: parameter a of the Angstrom model (intercept)
        angstr_b: parameter b of the Angstrom model (slope)
        
    Returns:
        Disaggregated hourly values of shortwave radiation.
    """
    # check if disaggregation method has a valid value
    if method not in ('pot_rad', 'pot_rad_via_ssd', 'pot_rad_via_bc'):
        raise ValueError('Invalid option')

    glob_disagg = pd.Series(index=melodist.util.hourly_index(data_daily.index))
    pot_rad_daily = pot_rad.resample('D').mean()

    if method == 'pot_rad':
        globalrad = data_daily.glob
    elif method == 'pot_rad_via_ssd':
        # in this case use the Angstrom model
        globalrad = pd.Series(index=data_daily.index, data=0.)
        dates = sun_times.index[sun_times.daylength > 0] # account for polar nights
        globalrad[dates] = angstroem(data_daily.ssd[dates], sun_times.daylength[dates], pot_rad_daily[dates], angstr_a, angstr_b)
    elif method == 'pot_rad_via_bc':
        # using data from Bristow-Campbell model
        globalrad = bristow_campbell(data_daily.tmin, data_daily.tmax, pot_rad_daily, bristcamp_a, bristcamp_c)

    globalrad_equal = globalrad.reindex(pot_rad.index, method='ffill') # hourly values (replicate daily mean value for each hour)
    pot_rad_daily_equal = pot_rad_daily.reindex(pot_rad.index, method='ffill')
    glob_disagg = pot_rad / pot_rad_daily_equal * globalrad_equal
    glob_disagg[glob_disagg < 1e-2] = 0.

    return glob_disagg


def potential_radiation(dates, lon, lat, timezone):    
    """calculates potential shortwave radiation for apecific location and time

    This routine calculates global radiation as described in:
    Liston, G. E. and Elder, K. (2006): A Meteorological Distribution System for High-Resolution Terrestrial Modeling (MicroMet), 
    J. Hydrometeorol., 7, 217–234.
    
    Corrections for eccentricity are carried out following:
    Paltridge, G.W., Platt, C.M.R., 1976. Radiative processes in Meteorology and Climatology. 
    Elsevier Scientific Publishing Company, Amsterdam, Oxford, New York.

    Args:
        dates: series (index) for which potential radiation shall be calculated
        lon: longitude (degree)
        lat: latitude (degree)
        timezone: timezone shift between local time and UTC
    Returns:
        series of potential shortwave radiation
    """
    

    solar_constant = 1367.
    days_per_year   = 365.25
    tropic_of_cancer = 0.41
    solstice = 173.0
    cloud_fraction = 0.

    glob = pd.Series(index=dates)

    for datetime in dates:
        day_of_year = (datetime - pd.datetime(datetime.year, 1, 1)).days + 1

        # compute solar decline in rad
        solar_decline = tropic_of_cancer * cos(2.0 * pi * (float(day_of_year) - solstice) / days_per_year)

        # compute the sun hour angle in rad
        standard_meridian = timezone * 15.
        delta_lat_time = (lon - standard_meridian) * 24. / 360.
        hour_angle = np.pi * (((datetime.hour + datetime.minute / 60. + delta_lat_time) / 12.) - 1.)

        # get solar zenit angle 
        cos_solar_zenith = sin(solar_decline) * sin(np.deg2rad(lat)) + cos(solar_decline) * cos(np.deg2rad(lat)) * cos(hour_angle)

        if cos_solar_zenith < 0:
            cos_solar_zenith = 0.0

        solar_zenith_angle = acos(cos_solar_zenith)

        # compute transmissivities for direct and diffus radiation using cloud fraction
        transmissivity_direct = (0.6 + 0.2 * cos_solar_zenith) * (1.0 - cloud_fraction)
        transmissivity_diffuse = (0.3 + 0.1 * cos_solar_zenith) * cloud_fraction

        # modify solar constant for eccentricity
        beta = 2. * pi * (day_of_year / days_per_year)
        radius_ratio = 1.00011 + 0.034221 * cos(beta) + 0.00128 * sin(beta) + 0.000719 * cos(2. * beta) + 0.000077 * sin(2 * beta)
        solar_constant_times_radius_ratio = solar_constant * radius_ratio

        # get total shortwave radiation
        direct_radiation = (solar_constant_times_radius_ratio * transmissivity_direct) * cos(solar_zenith_angle)
        diffuse_radiation = solar_constant_times_radius_ratio * transmissivity_diffuse # * cos(solarZenitAngle) #?

        if direct_radiation < 0:
            diffuse_radiation = 0
            direct_radiation = 0.

        glob[datetime] = direct_radiation + diffuse_radiation

    return glob


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
    temp = temp.reindex(pd.DatetimeIndex(start=temp.index[0], end=temp.index[-1], freq='D'))
    temp['tmin_nextday'] = temp.tmin
    temp.tmin_nextday.iloc[:-1] = temp.tmin.iloc[1:].values

    temp = temp.loc[tmin.index]
    pot_rad_daily = pot_rad_daily.loc[tmin.index]

    dT = temp.tmax - (temp.tmin + temp.tmin_nextday) / 2

    dT_m_avg = dT.groupby(dT.index.month).mean()
    B = 0.036 * np.exp(-0.154 * dT_m_avg[temp.index.month])
    B.index = temp.index

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
    df = pd.DataFrame(data=dict(tmin=tmin, tmax=tmax, pot=pot_rad_daily, obs=obs_rad_daily)).dropna(how='any')
    bc_absbias = lambda ac: np.abs(np.mean(bristow_campbell(df.tmin, df.tmax, df.pot, ac[0], ac[1]) - df.obs))
    res = scipy.optimize.minimize(bc_absbias, [0.75, 2.4]) # i.e. we minimize the absolute bias

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

    angstroem_opt = lambda x, a, b: angstroem(x[0], x[1], x[2], a, b)

    x = np.array([df.ssd, df.day_length, df.pot])
    popt, pcov = scipy.optimize.curve_fit(angstroem_opt, x, df.obs, p0=[0.25, 0.75])

    return popt
