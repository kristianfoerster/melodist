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

import pandas as pd
from datetime import timedelta
import glob


def get_variables():
    """returns the data fields relevant for MELODIST

    Returns:
        Array including names of meteorological variables
    """
    [
        'temp',
        'precip',
        'glob',
        'hum',
        'wind',
        'ssd',
    ]


def read_single_knmi_file(filename):
    """reads a single file of KNMI's meteorological time series

    data availability: www.knmi.nl/nederland-nu/klimatologie/uurgegevens

    Args:
        filename: the file to be opened

    Returns:
        pandas data frame including time series
    """
    hourly_data_obs_raw = pd.read_csv(
        filename,
        parse_dates=[['YYYYMMDD', 'HH']],
        date_parser=lambda yyyymmdd, hh: pd.datetime(int(str(yyyymmdd)[0:4]),
                                                     int(str(yyyymmdd)[4:6]),
                                                     int(str(yyyymmdd)[6:8]),
                                                     int(hh) - 1),
        skiprows=31,
        skipinitialspace=True,
        na_values='',
        keep_date_col=True,
    )

    hourly_data_obs_raw.index = hourly_data_obs_raw['YYYYMMDD_HH']
    hourly_data_obs_raw.index = hourly_data_obs_raw.index + timedelta(hours=1)

    columns_hourly = get_variables()

    hourly_data_obs = pd.DataFrame(
        index=hourly_data_obs_raw.index,
        columns=columns_hourly,
        data=dict(
            temp=hourly_data_obs_raw['T'] / 10 + 273.15,
            precip=hourly_data_obs_raw['RH'] / 10,
            glob=hourly_data_obs_raw['Q'] * 10000 / 3600.,
            hum=hourly_data_obs_raw['U'],
            wind=hourly_data_obs_raw['FH'] / 10,
            ssd=hourly_data_obs_raw['SQ'] * 6,
        ),
    )
    # remove negative values
    negative_values = hourly_data_obs['precip'] < 0.0
    hourly_data_obs.loc[negative_values, 'precip'] = 0.0
    return hourly_data_obs


def read_knmi_dataset(directory):
    """Reads files from a directory and merges the time series

    Please note: For each station, a separate directory must be provided!
    data availability: www.knmi.nl/nederland-nu/klimatologie/uurgegevens

    Args:
        directory: directory including the files

    Returns:
        pandas data frame including time series
    """
    filemask = '%s*.txt' % directory
    filelist = glob.glob(filemask)

    columns_hourly = get_variables()
    ts = pd.DataFrame(columns=columns_hourly)

    first_call = True
    for file_i in filelist:
        print(file_i)
        current = read_single_knmi_file(file_i)
        if(first_call):
            ts = current
            first_call = False
        else:
            ts = pd.concat([ts, current])
    return ts
