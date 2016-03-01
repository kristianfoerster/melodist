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
from math import cos, sin, acos, pi

"""
This routine diaggregates daily values of global radiation data to hourly values
"""


def disaggregate_radiation(data_daily, sun_times, pot_rad, method='pot_rad', angstr_a=0.25, angstr_b=0.5):
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
    pot_rad_daily = pot_rad.resample('D', how='sum')

    # call Bristow-Campbell model prior to radiation computations if required
    if method == 'pot_rad_via_bc':
        rad_bc = bristow_campbell(data_daily, pot_rad_daily)

    # for this option calculate incoming solar radiation as a function of the station location, time and cloud cover
    for index, row in data_daily.iterrows():
        # now account for the available radiation from the daily observations
        if method == 'pot_rad':
            globalrad = data_daily.glob[index]
        elif method == 'pot_rad_via_ssd':
            # in this case use the Angstrom model...
            # first step: calculate maximum sunshine duration

            # @todo: add time_zone to global radiation calculation as well
            # Call Angstrom model using the daily potential radiation as computed previously
            if sun_times.daylength[index] > 0:
                globalrad = (angstr_a + angstr_b * data_daily.ssd[index] / sun_times.daylength[index]) * pot_rad_daily[index] / 24
            else:
                globalrad = 0 # polar night case
        elif method == 'pot_rad_via_bc':
            # using data from Bristow-Campbell model
            globalrad = rad_bc[index] / 24

        for hour in range(0, 24):
            datetime = index.replace(hour=hour)

            # error handling: if globalrad is invalid, apply 0.5 * potential value as a crude assumption
            # if np.isnan(globalrad):
            #     globalrad = 0.5 * pot_rad_daily[index] / 24
            if pot_rad_daily[index] > 0.:
                glob_disagg[datetime] = (pot_rad[datetime]/pot_rad_daily[index]) * (globalrad * 24.)
            else:
                glob_disagg[datetime] = 0.  # handle polar night values

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


def bristow_campbell(data_daily, pot_rad, A = 0.75, C = 2.4):
    """calculates potential shortwave radiation based on minimum and maximum temperature

    This routine calculates global radiation as described in:
    Bristow, Keith L., and Gaylon S. Campbell: On the relationship between
    incoming solar radiation and daily maximum and minimum temperature.
    Agricultural and forest meteorology 31.2 (1984): 159-166.

    Args:
        daily_data: time series (daily data) including at least minimum and maximum temeprature
        pot_rad: hourly dataframe including potential radiation
        A: parameter A of the Bristow-Campbell model
        C: parameter C of the Bristow-Campbell model
    Returns:
        series of potential shortwave radiation
    """


    # nmean indicates the number of days for calculating the average difference between Tn and Tx    
    R0 = pd.Series(index=data_daily.index)
    n_steps = len(R0)
    avDeltaT = np.zeros(n_steps)
    transmissivity = np.zeros(n_steps)

    # calculate mean daily temperature range for each month
    dT = data_daily.tmax - data_daily.tmin
    dT_monthly = dT.resample(rule="24h", how="mean")
    dT_m_avg = dT_monthly.groupby(dT_monthly.index.month).aggregate("mean")

    dT[np.isnan(dT)] = 0
    for i in range(0, n_steps):
        if i < (n_steps - 1) and ~np.isnan(data_daily.tmax[i]) and ~np.isnan(data_daily.tmin[i]) and ~np.isnan(data_daily.tmin[i+1]):
            avDeltaT[i] = data_daily.tmax[i] - 0.5 * (data_daily.tmin[i] + data_daily.tmin[i + 1]) # Tx(i) - 0.5* (Tn(i)+Tn(i+1));
        else:
            avDeltaT[i]=0.

        B = 0.036 * np.exp(-0.154 * dT_m_avg[data_daily.index.month[i]])
        transmissivity[i] = A * (1 - np.exp(-B * dT[i]**C))

        R0[i] = transmissivity[i] * pot_rad[i]
    return R0
