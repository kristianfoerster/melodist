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
from . import cascade
from .precipitation import *
from .temperature import *
from .radiation import *
from .humidity import *
from .wind import *
from .data_reader import *
from .util.util import distribute_equally
from .station import Station
from .stationstatistics import StationStatistics

columns_daily=[
    'tmean',
    'tmin',
    'tmax',
    'precip',
    'glob',
    'ssd',
    'hum',
    'wind',
]

columns_hourly=[
    'temp',
    'precip',
    'glob',
    'hum',
    'wind',
]

class Options: # (deprecated)
    TEMP_SINE_CURVE = 'sine'
    TEMP_LINEAR = 'linear'
    EQUAL_DISTRIBUTION = 'equal'
    PREC_CASCADE = 'cascade'
    PREC_MASTER_STATION = 'masterstation'
    HUMIDITY_MINIMAL = 'minimal'
    HUMIDITY_DEWPOINT_REGRESSION = 'dewpoint_regression'
    HUMIDITY_MIN_MAX = 'min_max'
    HUMIDITY_LINEAR_DEWPOINT_VARIATION = 'linear_dewpoint_variation'
    WIND_FITTED_COSINE = 'cosine'
    WIND_RANDOM_EXPONENTIAL = 'random'
    POT_RAD = 'pot_rad'
    POT_RAD_VIA_SSD = 'pot_rad_via_ssd'
    POT_RAD_VIA_BC = 'pot_rad_via_bc'
