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
import read_knmi_data

# This an example that demonstrates the usage of MELODIST.

# Step 1: load some data for testing:
#
# The data provided by the KNMI are well suited for testing. You might test
# MELODIST using the data for De Bilt.
# The entire catalog can be browsed for more data:
# https://www.knmi.nl/nederland-nu/klimatologie/uurgegevens
#
# download and read files for De Bilt (ID=260)
path_inp = 'https://cdn.knmi.nl/knmi/map/page/klimatologie/gegevens/uurgegevens/uurgeg_260_2001-2010.zip'

data_obs_hourly = read_knmi_data.read_single_knmi_file(path_inp)

# recommendation: shift datetime to local time zone
data_obs_hourly = data_obs_hourly.shift(1)  # UTC -> CET

# truncate dataframe in order to obtain full days if applicable
data_obs_hourly = melodist.util.drop_incomplete_days(data_obs_hourly)

# optional: use a subselection for calculating the various statistics
data_obs_hourly_sel = data_obs_hourly.loc['2009-01-01':'2009-12-31']

# for testing purposes in this file: create daily out of hourly values for
# disaggregation
data_obs_daily = melodist.util.daily_from_hourly(data_obs_hourly_sel)

# Create station object
longitude = 5.18
latitude = 52.10
timezone = 1
debilt = melodist.Station(lon=longitude, lat=latitude, timezone=timezone)
debilt.data_daily = data_obs_daily

# initialize station statistics object
debilt.statistics = melodist.StationStatistics(data_obs_hourly_sel)

# Step 2: Compute statistics for various methods
stats = debilt.statistics
stats.calc_wind_stats()
stats.calc_humidity_stats()
stats.calc_temperature_stats()
stats.calc_radiation_stats()

# the computations of statistics is similar:
stats.calc_precipitation_stats()
# as an alternative you can also define 2 seasons for precipitation statistics
# months_season1 = np.array([11,12,1,2,3,4]) # winter season
# months_season2 = np.array([5,6,7,8,9,10]) # summer season
# seasons = [months_season1, months_season2] # combine subsets to one array
# stats.calc_precipitation_stats(months=seasons)

# Step 3: disaggregate temperature
debilt.disaggregate_temperature(method='sine_min_max', min_max_time='fix')
# temp_sim_T1a = debilt.data_disagg.temp.copy()
# s.disaggregate_temperature(method='sine_min_max', min_max_time='sun_loc')
# temp_sim_T1b = debilt.data_disagg.temp.copy()
# s.disaggregate_temperature(method='sine_min_max', min_max_time='sun_loc_shift')
# temp_sim_T1c = debilt.data_disagg.temp.copy()
# s.disaggregate_temperature(method='sine_min_max', min_max_time='sun_loc', mod_nighttime=True)
# temp_sim_T1d = debilt.data_disagg.temp.copy()
# s.disaggregate_temperature(method='sine_min_max', min_max_time='sun_loc')

# Step 4: disaggregate humidity
debilt.disaggregate_humidity(method='linear_dewpoint_variation')
# hum_sim_H3 = debilt.data_disagg.hum.copy()
# debilt.disaggregate_humidity(method='dewpoint_regression')
# hum_sim_H2 = debilt.data_disagg.hum.copy()
# debilt.disaggregate_humidity(method='minimal')
# hum_sim_H1 = debilt.data_disagg.hum.copy()
# debilt.disaggregate_humidity(method='min_max')
# hum_sim_H4 = debilt.data_disagg.hum.copy()

# Step 5: disaggregate wind speed
debilt.disaggregate_wind(method='cosine')

# Step 6: disaggregate radiation
# Suppose we want to apply the Angstrom model to compute shortwave radiation.
# In this case, the Angstrom parameters have been calculated from the hourly data using
# stats.calc_radiation_stats(). They could also be supplied manually:
# debilt.statistics.glob.angstroem_a = 0.36  # valid for De Bilt
# debilt.statistics.glob.angstroem_b = 0.68  # valid for De Bilt
debilt.disaggregate_radiation(method='pot_rad_via_ssd')

# Step 7: disaggregate precipitation
debilt.disaggregate_precipitation(method='cascade')

# Step 8: Testing the methods for longer time series
# We have calculated the statistics from just one year of hourly data. Now we can read in
# the daily data and disaggregate it using the statistics stored in the Station object.
data_obs_daily = melodist.util.daily_from_hourly(data_obs_hourly)
# in this example we read the hourly data (all 10 years now) and manually resample it
# to daily values. Usually obviously we would read in actual daily values.
debilt.data_daily = data_obs_daily

# Now we can use the usual disaggregation routines, like in the examples before:
# Warning: the dataframes will be extended according to the new length of input data!
debilt.disaggregate_precipitation(method='cascade')

# Alternatively, we could also export the calculated statistics to a JSON file,
# create a new station, and import the statistics file:
debilt.statistics.to_json('/tmp/debilt_stats.json')
debilt2 = melodist.Station(lon=longitude, lat=latitude, timezone=timezone)
debilt2.data_daily = data_obs_daily
debilt2.statistics = melodist.StationStatistics.from_json('/tmp/debilt_stats.json')
