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
import numpy as np
import pandas as pd
import scipy.stats

def hourly_index(daily_index, fill_gaps=False):
    index = pd.DatetimeIndex(start=daily_index.min(), end=daily_index.max().replace(hour=23), freq='H')

    # remove days that are not in the daily index:
    if not fill_gaps and len(index) > len(daily_index) * 24:
        df = pd.DataFrame(index=index.date, data=dict(hour=index.hour, keep=True))
        df.index = df.index.to_datetime()
        dropdates = list(set(index.date) - set(daily_index.date))
        df.loc[dropdates, 'keep'] = False
        df = df[df.keep]
        index = pd.DatetimeIndex([pd.datetime(y, m, d, h) for y, m, d, h in zip(df.index.year, df.index.month, df.index.day, df.hour)])

    return index

def distribute_equally(daily_data, divide=False):
    """Obtains hourly values by equally distributing the daily values.

    Args:
        daily_data: daily values
        divide: if True, divide resulting values by the number of hours in
            order to preserve the daily sum (required e.g. for precipitation).

    Returns:
        Equally distributed hourly values.
    """

    index = hourly_index(daily_data.index)
    hourly_data = daily_data.reindex(index)
    hourly_data = hourly_data.groupby(hourly_data.index.day).transform(
        lambda x: x.fillna(method='ffill'))

    if divide:
        hourly_data /= 24

    return hourly_data

def vapor_pressure(temp, hum):
    """
    Calculates vapor pressure from temperature and humidity after Sonntag (1990).

    Args:
        temp: temperature values
        hum: humidity value(s). Can be scalar (e.g. for calculating saturation vapor pressure).

    Returns:
        Vapor pressure in hPa.
    """

    if np.isscalar(hum):
        hum = np.zeros(temp.shape) + hum

    assert(temp.shape == hum.shape)

    positives = np.array(temp >= 273.15)
    vap_press = np.zeros(temp.shape) * np.nan
    vap_press[positives] = 6.112 * np.exp((17.62 * (temp[positives] - 273.15)) / (243.12 + (temp[positives] - 273.15))) * hum[positives] / 100.
    vap_press[~positives] = 6.112 * np.exp((22.46 * (temp[~positives] - 273.15)) / (272.62  + (temp[~positives] - 273.15))) * hum[~positives] / 100.

    return vap_press

def dewpoint_temperature(temp, hum):
    """computes the dewpoint temperature

    Parameters
    ----
    temp :      temperature [K]
    hum :       relative humidity
    

    Returns
        dewpoint temperature in K
    """
    assert(temp.shape == hum.shape)

    vap_press = vapor_pressure(temp, hum)

    positives = np.array(temp >= 273.15)
    dewpoint_temp = temp.copy() * np.nan
    dewpoint_temp[positives] = 243.12 * np.log(vap_press[positives] / 6.112) / (17.62 - np.log(vap_press[positives] / 6.112))
    dewpoint_temp[~positives] = 272.62 * np.log(vap_press[~positives] / 6.112) / (22.46 - np.log(vap_press[~positives] / 6.112))

    return dewpoint_temp + 273.15

def linregress(x, y, return_stats=False):
    """linear regression calculation

    Parameters
    ----
    x :         independent variable (series)
    y :         dependent variable (series)
    return_stats : returns statistical values as well if required (bool)
    

    Returns
    ----
    list of parameters (and statistics)
    """
    a1, a0, r_value, p_value, stderr = scipy.stats.linregress(x, y)

    retval = a1, a0
    if return_stats:
        retval += r_value, p_value, stderr

    return retval

def get_sun_times(dates, lon, lat, time_zone):
    """Computes the times of sunrise, solar noon, and sunset for each day.

    Parameters
    ----
    dates:      datetime
    lat :       latitude in DecDeg
    lon :       longitude in DecDeg
    time_zone : timezone
    

    Returns
    ----
    DataFrame:  [sunrise, sunnoon, sunset, day length] in dec hours
    """

    df = pd.DataFrame(index=dates, columns=['sunrise', 'sunnoon', 'sunset', 'daylength'])

    doy = np.array([(d - d.replace(day=1, month=1)).days + 1 for d in df.index]) # day of year

    # Day angle and declination after Bourges (1985):
    day_angle_b = np.deg2rad((360. / 365.25) * (doy - 79.346))
    
    declination = np.deg2rad(
        0.3723 + 23.2567 * np.sin(  day_angle_b) - 0.7580 * np.cos(  day_angle_b)
            +  0.1149 * np.sin(2*day_angle_b) + 0.3656 * np.cos(2*day_angle_b)
            -  0.1712 * np.sin(3*day_angle_b) + 0.0201 * np.cos(3*day_angle_b)
            )
    
    # Equation of time with day angle after Spencer (1971):
    day_angle_s = 2 * np.pi * (doy - 1) / 365.
    eq_time = 12. / np.pi * (
        0.000075 +
        0.001868 * np.cos(  day_angle_s) - 0.032077 * np.sin(  day_angle_s) -
        0.014615 * np.cos(2*day_angle_s) - 0.040849 * np.sin(2*day_angle_s)
        )
    
    #
    standard_meridian = time_zone * 15.
    delta_lat_time = (lon - standard_meridian) * 24. / 360.
    
    omega_nul_arg = -np.tan(np.deg2rad(lat)) * np.tan(declination)
    omega_nul = np.arccos(omega_nul_arg)
    sunrise  = 12. * (1. - (omega_nul) / np.pi) - delta_lat_time - eq_time
    sunset   = 12. * (1. + (omega_nul) / np.pi) - delta_lat_time - eq_time

    # as an approximation, solar noon is independent of the below mentioned
    # cases:
    sunnoon  = 12. * (1.) - delta_lat_time - eq_time        
    
    # $kf 2015-11-13: special case midnight sun and polar night
    # CASE 1: MIDNIGHT SUN
    # set sunrise and sunset to values that would yield the maximum day
    # length even though this a crude assumption
    pos = omega_nul_arg < -1
    sunrise[pos] = sunnoon[pos] - 12
    sunset[pos]  = sunnoon[pos] + 12

    # CASE 2: POLAR NIGHT
    # set sunrise and sunset to values that would yield the minmum day
    # length even though this a crude assumption
    pos = omega_nul_arg > 1
    sunrise[pos] = sunnoon[pos]
    sunset[pos]  = sunnoon[pos]

    daylength = sunset - sunrise
        
    # adjust if required
    sunrise[sunrise < 0] += 24
    sunset[sunset > 24] -= 24

    df.sunrise = sunrise
    df.sunnoon = sunnoon
    df.sunset = sunset
    df.daylength = daylength

    return df

def detect_gaps(dataframe, timestep, print_all=False, print_max=5, verbose=True):
    """checks if a given dataframe contains gaps and returns the number of gaps

    This funtion checks if a dataframe contains any gaps for a given temporal
    resolution that needs to be specified in seconds. The number of gaps
    detected in the dataframe is returned.

    Args:
        dataframe: A pandas dataframe object with index defined as datetime
        timestep (int): The temporal resolution of the time series in seconds
            (e.g., 86400 for daily values)
        print_all (bool, opt): Lists every gap on the screen
        print_mx (int, opt): The maximum number of gaps listed on the screen in
            order to avoid a decrease in performance if numerous gaps occur
        verbose (bool, opt): Enables/disables output to the screen

    Returns:
        The number of gaps as integer. Negative values indicate errors.
    """
    gcount = 0
    msg_counter = 0
    warning_printed = False
    try:
        n = len(dataframe.index)
    except:
        print('Error: Invalid dataframe.')
        return -1
    for i in range(0,n):
        if(i > 0):
            time_diff = dataframe.index[i] - dataframe.index[i-1]
            if(time_diff.delta/1E9 != timestep):
                gcount += 1
                if(print_all or (msg_counter <= print_max - 1)):
                    if verbose:
                        print('Warning: Gap in time series found between %s and %s' % (dataframe.index[i-1], dataframe.index[i]))
                    msg_counter += 1
                if(msg_counter == print_max and warning_printed == False and verbose):
                    print('Waring: Only the first %i gaps have been listed. Try to increase print_max parameter to show more details.' % msg_counter)
                    warning_printed = True
    if verbose:
        print('%i gaps found in total.' % (gcount))
    return gcount

def drop_incomplete_days(dataframe, shift=0):
    """truncates a given dataframe to full days only

    This funtion truncates a given pandas dataframe (time series) to full days
    only, thus dropping leading and tailing hours of incomplete days. Please
    note that this methodology only applies to hourly time series.

    Args:
        dataframe: A pandas dataframe object with index defined as datetime
        shift (unsigned int, opt): First hour of daily recordings. For daily
            recordings of precipitation gages, 8 would be the first hour of
            the subsequent day of recordings since daily totals are
            usually recorded at 7. Omit defining this parameter if you intend
            to pertain recordings to 0-23h.

    Returns:
        A dataframe with full days only.
    """
    dropped = 0
    if shift > 23 or shift < 0:
        print("Invalid shift parameter setting! Using defaults.")
        shift = 0
    first = shift
    last  = first - 1
    if last < 0:
        last += 24
    try:
        # todo: move this checks to a separate function
        n = len(dataframe.index)
    except:
        print('Error: Invalid dataframe.')
        return dataframe
    
    delete = list()  
    
    # drop heading lines if required
    for i in range(0,n):
        if dataframe.index.hour[i] == first and dataframe.index.minute[i] == 0:
            break
        else:
            delete.append(i)
            dropped += 1

    # drop tailing lines if required
    for i in range(n-1,0,-1):
        if dataframe.index.hour[i] == last and dataframe.index.minute[i] == 0:
            break
        else:
            delete.append(i)
            dropped += 1
    # print("The following rows have been dropped (%i in total):" % dropped)
    # print(delete)
    return dataframe.drop(dataframe.index[[delete]])

def prepare_interpolation_data(data_daily, column_hours):
    start_date = data_daily.index[0]
    end_date = data_daily.index[-1]

    data = pd.Series()

    for column, hour in column_hours.items():
        index = pd.DatetimeIndex(start=start_date.replace(hour=hour), end=end_date.replace(hour=hour), freq='D')
        s = pd.Series(index=index, data=data_daily[column].values)
        data = data.append(s)

    data = data.sort_index()
    data = data.reindex(hourly_index(data_daily.index))

    return data

def daily_from_hourly(df):
    """Aggregates data (hourly to daily values) according to the characteristics
    of each variable (e.g., average for temperature, sum for precipitation)

    Args:
        df: dataframe including time series with one hour time steps

    Returns:
        dataframe (daily)
    """

    df_daily = pd.DataFrame()

    if 'temp' in df:
        df_daily['temp'] = df.temp.resample('D').mean()
        df_daily['tmin'] = df.temp.groupby(df.temp.index.date).min()
        df_daily['tmax'] = df.temp.groupby(df.temp.index.date).max()

    if 'precip' in df:
        df_daily['precip'] = df.precip.resample('D').sum()

    if 'glob' in df:
        df_daily['glob'] = df.glob.resample('D').mean()

    if 'hum' in df:
        df_daily['hum'] = df.hum.resample('D').mean()

    if 'hum' in df:
        df_daily['hum_min'] = df.hum.groupby(df.hum.index.date).min()

    if 'hum' in df:
        df_daily['hum_max'] = df.hum.groupby(df.hum.index.date).max()

    if 'wind' in df:
        df_daily['wind'] = df.wind.resample('D').mean()

    if 'ssd' in df:
        df_daily['ssd'] = df.ssd.resample('D').sum() / 60 # minutes to hours

    df_daily.index.name = None
    return df_daily

