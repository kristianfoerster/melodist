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
from . __version__ import __version__

from . import cascade
from .data_io import *
from .humidity import *
from .precipitation import *
from .radiation import *
from .station import Station
from .stationstatistics import StationStatistics
from .temperature import *
from .util.util import distribute_equally
from .wind import *

columns_daily = [
    'tmean',
    'tmin',
    'tmax',
    'precip',
    'glob',
    'ssd',
    'hum',
    'wind',
]

columns_hourly = [
    'temp',
    'precip',
    'glob',
    'hum',
    'wind',
]


class Options:  # (deprecated)
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
