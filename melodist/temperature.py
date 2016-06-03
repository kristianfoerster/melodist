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
import melodist.util
import numpy as np

"""
This routine diaggregates daily values of temperature data to hourly values
"""

import re
import pandas as pd
from math import sin, cos, pi


def disaggregate_temperature(dailyData, method='sine', min_max_time='fix', mod_nighttime = False, max_delta=None, sun_times=None):
    """The disaggregation function for temperature

    Parameters
    ----
    dailyData :       daily data
    method :          method to dissaggregate
    min_max_time:     "fix" - min/max temperature at fix times 7h/14h,
                      "sun_loc" - min/max calculated by sunrise/sunnoon + 2h,
                      "sun_loc_shift" - min/max calculated by sunrise/sunnoon + monthly mean shift ,                     
    max_delta:        maximum monthly temperature shift as returned by get_shift_by_data()
    sun_times:        times of sunrise/noon as returned by get_sun_times()
    """

    # get length of time series in days
    days = len(dailyData.index)

    # perform disaggregation depending on chosen option

    if method not in (
        'sine',
    ):
        raise ValueError('Invalid option')

    temp_disagg = pd.Series(index=melodist.util.hourly_index(dailyData.index))

###############################################################################

    # for this option assume maximum occur at 14h and minimm at 7h and fit cosine function through minimum and maximum temperatures
    if method == 'sine':
        # set some parameters for cosine function
        steps = 24

        # init first value for polar night disaggregation
        last_val = 0.5 * (dailyData.tmin[[0]] + dailyData.tmin[[0]])
        model_transition = False
        daylength_thres = 3
        # min / max hour during polar night assumption
        minloc = 6
        maxloc = 18
 

        # loop over all days in time series
        counter = 0
        for index, row in dailyData.iterrows():
            
            #take fix location for minimum and maximum
            if min_max_time== "fix":
                minLocation= 7
                maxLocation= 14
        
            #take location for minimum and maximum by sunrise / sunnoon + 2h
            elif min_max_time== "sun_loc":
                minLocation= int(round(sun_times.sunrise[index]))      #sun rise round to full hour 
                maxLocation= int(round(sun_times.sunnoon[index])) + 2  #sun noon round to full hour + fix 2h
            
            #take location for minimum and maximum by sunrise / sunnoon + monthly delta
            elif min_max_time== "sun_loc_shift":
                minLocation= int(round(sun_times.sunrise[index]))                  #sun rise round to full hour 
                maxLocation= int(round(sun_times.sunnoon[index] + max_delta[index.month]))  #sun noon + shift derived from obsereved h data, round to full hour
                if minLocation > maxLocation:
                    maxLocation= int(round(sun_times.sunnoon[index])) + 2 # standard shift in this case
                # month -1 ???
                
            # define minimum value for cosine function to assure smooth transitions
            # for last timestep there is no "next step" so consider two cases in minValue calculation
            if counter == days - 1 :   
                minValueAct = row.tmin
                minValueNext = row.tmin
                maxValueNext = row.tmax
                
            else:
                minValueAct = row.tmin
                minValueNext = dailyData.tmin[index.date() + pd.Timedelta('1 days')] #get min value of next day
                maxValueNext = dailyData.tmax[index.date() + pd.Timedelta('1 days')]


            # define maximum value for cosine function to assure smooth transitions
            # for first timestep there is no "step before" so consider two cases in maxValue calculation
            if counter == 0 :
                maxValueAct = row.tmax
                maxValueBefore = row.tmax
                minValueBefore = row.tmin

            else:
                maxValueAct = row.tmax
                maxValueBefore = dailyData.tmax[index.date() - pd.Timedelta('1 days')] #get max value of day before
                minValueBefore = dailyData.tmin[index.date() - pd.Timedelta('1 days')] #get max value of day before
            
            # get value for current hour
            dayIndex = str(index)
            for hour in range(0, 24):
                if hour < 10:
                    hourIndex = dayIndex[0:11]+'0'+str(hour)+dayIndex[13:]
                else:
                    hourIndex = dayIndex[0:11]+str(hour)+dayIndex[13:]
    
                # during polar night, no diurnal variation of temperature is applied
                # instead the daily average calculated using tmin and tmax is applied
                if sun_times.daylength[index] < daylength_thres:
                    # temp_disagg[hourIndex] = 0.5 * (row.tmin + row.tmax)
                    avgBefore  = 0.5 * (minValueBefore + maxValueBefore)
                    avgCurrent = 0.5 * (minValueAct + maxValueAct)
                    avgNext    = 0.5 * (minValueNext + maxValueNext)
                    
                    maxloc_transition = maxloc                    
                    
                    if model_transition and hour <= minloc:
                        last_val = maxValueBefore
                        maxloc_transition = maxLocation
                        if hour == minloc:
                            model_transition = False
                                                                                        
                    if avgBefore <= avgCurrent:
                        val1 = minValueAct
                        val2 = maxValueAct
                    else:
                        val1 = maxValueAct
                        val2 = minValueAct

                    if hour <= minloc:
                        funcValue = last_val + (hour + 24 - maxloc_transition)/ (24 - maxloc_transition + minloc) * (val1 - last_val)
                    elif hour > minloc and hour <= maxloc:
                        funcValue = val1 + (hour - minloc)/ (maxloc - minloc) * (val2 - val1)
                    else:
                        if avgCurrent <= avgNext:
                            next_val = minValueNext
                        else:
                            next_val = maxValueNext

                        next_time = minloc
                        # check whether to use cosine function in next time step
                        if counter < days - 1:
                            if sun_times.daylength[counter+1] >= daylength_thres:
                                next_val = minValueNext
                                next_time = minLocation
                                if hour == 23:
                                    model_transition = True

                        funcValue = val2 + (hour - maxloc)/ (24 - maxloc + next_time) * (next_val - val2)
                        last_val = val2
                        

                else:
                    # whenever we are before the maximum for the current day, use minimum value of current day for cosine function fitting
                    # once we have passed the maximum value use the minimum for next day to ensure smooth transitions
                    if hour < maxLocation:
                        minValue = minValueAct
                    else:
                        minValue = minValueNext
    
                    # whenever we are before the minimum for the current day, use maximum value of day before for cosine function fitting
                    # once we have passed the minimum value use the maximum for the actual day to ensure smooth transitions
                    if hour < minLocation:
                        maxValue = maxValueBefore
                    else:
                        maxValue = maxValueAct               
                    
                    deltaValue = maxValue - minValue
                    vTrans = minValue + (0.5*deltaValue)
    
                    if mod_nighttime:
                        if hour <= minLocation:
                            funcValue = vTrans + (deltaValue*0.5) * cos(pi/(steps - (maxLocation-minLocation))*(steps - maxLocation + hour))
                        elif hour > minLocation and hour < maxLocation:
                            funcValue = vTrans + (deltaValue*0.5) * cos(1.25*pi + 0.75*pi/(maxLocation-minLocation)*(hour-minLocation))
                        else:
                            funcValue = vTrans + (deltaValue*0.5) * cos(pi/(steps - (maxLocation-minLocation))*(hour-maxLocation))
                    else:
                        funcValue = vTrans + (deltaValue*0.5) * cos(2*pi/steps*(hour-maxLocation))

                    # overwrite values for model transition
                    if model_transition and hour < minLocation:
                        funcValue = last_val + (hour + 24 - maxloc)/ (24 - maxloc + minLocation) * (minValueAct - last_val)

                    if model_transition and hour == minLocation:
                        model_transition = False
                    
                    # check whether to use cosine function in next time step
                    if counter < days - 1:
                        if sun_times.daylength[counter+1] < daylength_thres:
                            if hour == maxLocation:
                                avgCurrent = 0.5 * (minValueAct + maxValueAct)
                                avgNext    = 0.5 * (minValueNext + maxValueNext)
                                if avgCurrent < avgNext:
                                    val = minValueNext
                                else:
                                    val = maxValueNext
                            if hour > maxLocation:
                                funcValue = maxValueAct + (hour - maxLocation)/ (24 - maxLocation + minloc) * (val - maxValueAct)
                                if hour == 23:
                                    model_transition = True
                                    last_val = maxValueAct
                        
                    # correction when tmax current day < tmin next day
                    #if hour >= maxLocation: #if maxValueAct < minValueNext and hour >= maxLocation:
                    #    funcValue = maxValueAct + (hour - maxLocation) / (24 - maxLocation + minLocation) * (minValueNext - maxValueAct)

                    #if hour < minLocation: #if maxValueBefore < minValueAct and hour <= minLocation:
                    #    funcValue = maxValueBefore + (hour + 24 - maxLocation)/ (24 - maxLocation + minLocation) * (minValueAct - maxValueBefore)
                     
                temp_disagg[hourIndex] = funcValue

            counter = counter + 1

    return temp_disagg

###############################################################################
def get_shift_by_data(temp_hourly, lon, lat, time_zone):
    '''function to get max temp shift (monthly) by hourly data
    
    Parameters
    ----
    hourly_data_obs : observed hourly data 
    lat :             latitude in DezDeg
    lon :             longitude in DezDeg
    time_zone:        timezone
    '''
    #prepare a daily index 
    days = temp_hourly.resample('D').max()
    max_delta = days * np.nan

    sun_times = melodist.util.get_sun_times(days.index, lon, lat, time_zone)
    
    #get hourly data day by day
    for index_d, row in days.iteritems():        
        index = index_d.date().isoformat()
        temp_h = temp_hourly[index]
        
        if temp_h.empty or temp_h.isnull().any(): # hasnans:
            max_delta[index] = np.nan
        else:
            #get hour of max temp
            max_temp = temp_h.idxmax().hour
            
            #get sun min/max loction (dez. h)
            sun_maxLocation= (sun_times.sunnoon[index_d])
            
            #get delta
            delta_max = max_temp - sun_maxLocation
    
            #write to daily pd df
            max_delta[index] = delta_max
        
    months = max_delta.resample('M').mean()
    data_month_mean = months.groupby(months.index.month).agg('mean')
    shift_max_month_mean = data_month_mean.transpose()
    
    return shift_max_month_mean #max_delta
