# -*- coding: utf-8 -*-
########################################################################
# This file is part of MELODIST - MEteoroLOgical observation time      #
# series DISaggregation Tool a program to disaggregate daily values    #
# of meteorological variables to hourly values                         #
#                                                                      #
#                                                                      #
# Copyright (C) 2016  Florian Hanzer (1, 2), Kristian Förster (1, 2),  #
# Benjamin Winter (1, 2), Thomas Marke (1)                             #
#                                                                      #
# (1) Institute of Geography, University of Innsbruck, Austria         #
# (2) alpS - Centre for Climate Change Adaptation, Innsbruck, Austria  #
#                                                                      #
# MELODIST is free software: you can redistribute it and/or modify     #
# it under the terms of the GNU General Public License as published by #
# the Free Software Foundation, either version 3 of the License, or    #
# (at your option) any later version.                                  #
#                                                                      #
# MELODIST is distributed in the hope that it will be useful,          #
# but WITHOUT ANY WARRANTY; without even the implied warranty of       #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the         #
# GNU General Public License for more details.                         #
#                                                                      #
# You should have received a copy of the GNU General Public License    #
# along with this program.  If not, see <http://www.gnu.org/licenses/>.#
#                                                                      #
########################################################################

# Update 2018: Allow cascade-based disaggregation from daily data
# to 5 minutes data (9, 10, or 11 levels of branching)
# Author: Siling Chen

from __future__ import print_function, division, absolute_import

import numpy as np
import pandas as pd

import copy
import melodist
import melodist.util

from . import cascade


def disagg_prec(dailyData,
                method='equal',
                cascade_options=None,
                hourly_data_obs=None,
                zerodiv="uniform",
                shift=0):
    """The disaggregation function for precipitation.

    Parameters
    ----------
    dailyData : pd.Series
        daily data
    method : str
        method to disaggregate
    cascade_options : cascade object
        including statistical parameters for the cascade model
    hourly_data_obs : pd.Series
        observed hourly data of master station
    zerodiv : str
        method to deal with zero division by key "uniform" --> uniform
        distribution
    shift : int
        shifts the precipitation data by shift (int) steps (eg +7 for
        7:00 to 6:00)
    """

    if method not in ('equal', 'cascade', 'masterstation'):
        raise ValueError('Invalid option')

    if method == 'equal':
        precip_disagg = melodist.distribute_equally(dailyData.precip,
                                                    divide=True)
    elif method == 'masterstation':
        precip_disagg = precip_master_station(dailyData,
                                              hourly_data_obs,
                                              zerodiv)
    elif method == 'cascade':
        assert cascade_options is not None
        precip_disagg = disagg_prec_cascade(dailyData,
                                            cascade_options,
                                            shift=shift)

    return precip_disagg


def disagg_prec_cascade(precip_daily, 
                        cascade_options, hourly=True,level=9,
                        shift=0,
                        test=False):
    """Precipitation disaggregation with cascade model (Olsson, 1998)

    Parameters
    ----------
    precip_daily : pd.Series
        daily data
    hourly: Boolean (for an hourly resolution disaggregation)
        if False, then returns 5-min disaggregated precipitation 
        (disaggregation level depending on the "level" variable)
    cascade_options : cascade object
        including statistical parameters for the cascade model
    shift : int
        shifts the precipitation data by shift steps (eg +7 for 7:00 to
        6:00)
    test : bool
        test mode, returns time series of each cascade level
    """

    if len(precip_daily) < 2:
        raise ValueError('Input data must have at least two elements.')

    # set missing values to zero:
    precip_daily = precip_daily.copy()
    missing_days = precip_daily.index[precip_daily.isnull()]
    precip_daily[missing_days] = 0

    if hourly:
        si = 5  # index of first level
    else:
        si = level
    # statistics for branching into two bins
    wxxcum = np.zeros((7, 2, 4))

    if isinstance(cascade_options, melodist.cascade.CascadeStatistics):
        # this is the standard case considering one data set for all levels
        # get cumulative probabilities for branching
        overwrite_stats = False
        for k in range(0, 7):
            wxxcum[k, :, :] = cascade_options.wxx[k, :, :]
            if k > 0:
                wxxcum[k, :, :] = wxxcum[k-1, :, :] + wxxcum[k, :, :]
    elif isinstance(cascade_options, list):
        if len(cascade_options) == si:#5
            overwrite_stats = True
            list_casc = cascade_options
        else:
            raise ValueError('Cascade statistics list must have %s elements!' % si)
    else:
        raise TypeError('cascade_options has invalid type')

    # arrays for each level
    n = len(precip_daily)
    vdn1 = np.zeros(n*2)
    vdn2 = np.zeros(n*4)
    vdn3 = np.zeros(n*8)
    vdn4 = np.zeros(n*16)
    vdn5 = np.zeros(n*32)
    if not hourly:
        vdn6 = np.zeros(n*64)
        vdn7 = np.zeros(n*128)
        vdn8 = np.zeros(n*256)
        vdn9 = np.zeros(n*512)
    if level == 10 or level == 11:
        vdn10 = np.zeros(n*1024)
        if level == 11:
            vdn11 = np.zeros(n*2048)
            
    # class boundaries for histograms
    wclassbounds = np.array([0.0,
                             0.1429,
                             0.2857,
                             0.4286,
                             0.5714,
                             0.7143,
                             0.8571,
                             1.0])

    # disaggregation for each level
    for l in range(1, si+1):
        if l == 1:
            vdn_in = precip_daily
            vdn_out = vdn1
        elif l == 2:
            vdn_in = vdn_out
            vdn_out = vdn2
        elif l == 3:
            vdn_in = vdn_out
            vdn_out = vdn3
        elif l == 4:
            vdn_in = vdn_out
            vdn_out = vdn4
        elif l == 5:
            vdn_in = vdn_out
            vdn_out = vdn5
       
        elif l == 6:
            vdn_in = vdn_out
            vdn_out = vdn6
        elif l == 7:
            vdn_in = vdn_out
            vdn_out = vdn7
        elif l == 8:
            vdn_in = vdn_out
            vdn_out = vdn8
        elif l == 9:
            vdn_in = vdn_out
            vdn_out = vdn9
        elif l == 10:
            vdn_in = vdn_out
            vdn_out = vdn10
        elif l == 11:
            vdn_in = vdn_out
            vdn_out = vdn11
                
        si -= 1
        if overwrite_stats:
            cascade_options = list_casc[si]
            for k in range(0, 7):
                wxxcum[k, :, :] = cascade_options.wxx[k, :, :]
                if k > 0:
                    wxxcum[k, :, :] = wxxcum[k-1, :, :] + wxxcum[k, :, :]
            meanvol = cascade_options.threshold[0]
        else:
            meanvol = cascade_options.threshold[si]

        # evaluate mean rainfall intensity for wet boxes
        # these values should be determined during the aggregation phase!!!!!
        # mean volume threshold
        # meanvol = np.mean(vdn_in[vdn_in>0.])

        # use values derived parameter by parameter estimation instead
        # see above

        j = 0
        for i in range(0, len(vdn_in)):
            # it's raining now?
            if vdn_in[i] > 0:
                # determine type of box
                if i == 0:  # only starting or isolated
                    if vdn_in[i+1] > 0:
                        vbtype = cascade.BoxTypes.starting
                    else:
                        vbtype = cascade.BoxTypes.isolated
                elif i == len(vdn_in)-1:  # only ending or isolated
                    if vdn_in[i-1] > 0:
                        vbtype = cascade.BoxTypes.ending
                    else:
                        vbtype = cascade.BoxTypes.isolated
                else:  # neither at at the end nor at the beginning
                    if vdn_in[i-1] == 0 and vdn_in[i+1] == 0:
                        vbtype = cascade.BoxTypes.isolated

                    if vdn_in[i-1] == 0 and vdn_in[i+1] > 0:
                        vbtype = cascade.BoxTypes.starting

                    if vdn_in[i-1] > 0 and vdn_in[i+1] > 0:
                        vbtype = cascade.BoxTypes.enclosed

                    if vdn_in[i-1] > 0 and vdn_in[i+1] == 0:
                        vbtype = cascade.BoxTypes.ending

                # above or below mean?
                if vdn_in[i] > meanvol:
                    belowabove = 1  # above mean
                else:
                    belowabove = 0  # below mean

                #
                p = np.zeros((3, 1))
                p[0] = cascade_options.p01[belowabove, vbtype-1]  # index changed!
                p[1] = cascade_options.p10[belowabove, vbtype-1]
                p[2] = cascade_options.pxx[belowabove, vbtype-1]

                # draw a random number to determine the braching type
                rndp = np.random.random()

                if rndp <= p[0]:
                    # first box 0, second box: 1 P(0/1)
                    vdn_out[j] = 0.0
                    j = j + 1
                    vdn_out[j] = vdn_in[i]
                    j = j + 1
                elif rndp > p[0] and rndp <= p[0] + p[1]:
                    # first box 1, second box: 0 P(1/0)
                    vdn_out[j] = vdn_in[i]
                    j = j + 1
                    vdn_out[j] = 0.0
                    j = j + 1
                else:
                    # both boxes wet
                    # we need a new random number
                    rndw = np.random.random()

                    # guess w1:
                    for k in range(0, 7):
                        if rndw <= wxxcum[k, belowabove, vbtype-1]:
                            w1 = wclassbounds[k+1] - 1./14.  # class center
                            break

                    vdn_out[j] = w1 * vdn_in[i]
                    j = j + 1
                    vdn_out[j] = (1. - w1) * vdn_in[i]
                    j = j + 1

                    # check results (in the previous version this error has never been observed)
                    if w1 < 0 or w1 > 1:
                        print('error')
                        return
            else:
                # add two dry boxes
                vdn_out[j] = 0.0
                j = j + 1
                vdn_out[j] = 0.0
                j = j + 1
    
    if hourly:
        # uniformly disaggregate 0.75 h values to 0.25 h values
        vdn_025 = np.zeros(len(vdn_out)*3)
        j = 0
        for i in range(0, len(vdn_out)):
            for m in range(0, 3):
                vdn_025[j+m] = vdn_out[i] / 3.
            j = j + 3
    
        # aggregate to hourly time steps
        vdn_025cs = np.cumsum(vdn_025)
        vdn = np.zeros(int(len(vdn_025)/4))
        for i in range(0, len(vdn)+1):
            # for first hour take 4th item
            if i == 0:
                vdn[i] = vdn_025cs[3]
            elif i == 1:
                pass
            else:
                # >1 (starting with 2-1 = 1 item)
                vdn[i-1] = vdn_025cs[(i*4)-1] - vdn_025cs[(i*4)-5]
    
        disagg_precip = pd.Series(index=melodist.util.hourly_index(precip_daily.index), data=vdn)
    
    else:
        precip_sn = pd.Series(index= sub_level_index(precip_daily.index, level=level, fill_gaps=False), data=vdn_out)
        disagg_precip = precip_sn.resample('5min').sum()
    
    # set missing days to nan again:
    for date in missing_days:
         disagg_precip[ disagg_precip.index.date == date.date()] = np.nan

    # shifts the data by shift steps (fills with nan/cuts edge data )
    if shift != 0:
         disagg_precip =  disagg_precip.shift(shift) #? freq='1U')

    # return time series
    if test:
        if hourly:
            return vdn1, vdn2, vdn3, vdn4, vdn5, vdn_025, disagg_precip
        else:
            if level == 9:
                return vdn1, vdn2, vdn3, vdn4, vdn5, vdn6, vdn7, vdn8, vdn9, precip_sn, disagg_precip
            elif level == 10:
                return vdn1, vdn2, vdn3, vdn4, vdn5, vdn6, vdn7, vdn8, vdn9, vdn10, precip_sn, disagg_precip
            else:
                return vdn1, vdn2, vdn3, vdn4, vdn5, vdn6, vdn7, vdn8, vdn9, vdn10, vdn11, precip_sn, disagg_precip
    else:
        return  disagg_precip

def precip_master_station(precip_daily,
                          master_precip_hourly,
                          zerodiv):
    """Disaggregate precipitation based on the patterns of a master station

    Parameters
    -----------
    precip_daily : pd.Series
        daily data
    master_precip_hourly :  pd.Series
        observed hourly data of the master station
    zerodiv : str
        method to deal with zero division by key "uniform" --> uniform
        distribution
    """

    precip_hourly = pd.Series(index=melodist.util.hourly_index(precip_daily.index))

    # set some parameters for cosine function
    for index_d, precip in precip_daily.iteritems():

        # get hourly data of the day
        index = index_d.date().isoformat()
        precip_h = master_precip_hourly[index]

        # calc rel values and multiply by daily sums
        # check for zero division
        if precip_h.sum() != 0 and precip_h.sum() != np.isnan(precip_h.sum()):
            precip_h_rel = (precip_h / precip_h.sum()) * precip

        else:
            # uniform option will preserve daily data by uniform distr
            if zerodiv == 'uniform':
                precip_h_rel = (1/24) * precip

            else:
                precip_h_rel = 0

        # write the disaggregated day to data
        precip_hourly[index] = precip_h_rel

    return precip_hourly


def aggregate_precipitation(vec_data,hourly=True, percentile=50):
    """Aggregates highly resolved precipitation data and creates statistics

    Parameters
    ----------
    vec_data : pd.Series
        hourly (hourly=True) OR 5-min values 

    Returns
    -------
    output : cascade object
        representing statistics of the cascade model
    """
    cascade_opt = cascade.CascadeStatistics()
    cascade_opt.percentile = percentile

    # length of input time series
    n_in = len(vec_data)
    n_out = np.floor(n_in/2)

    # alternative:
    # 1st step: new time series
    vec_time = vec_data.index
    vdn0 = []
    vtn0 = []
    j = 0
    for i in range(0, n_in):
        if np.mod(i, 2) != 0:
            vdn0.append(vec_data.precip.values[i-1] + vec_data.precip.values[i])
            vtn0.append(vec_time[i])
            j = j+1

    vdn = pd.DataFrame(index=vtn0, data={'precip': vdn0})

    # length of new time series
    n_out = len(vdn)

    # series of box types:
    vbtype = np.zeros((n_out, ), dtype=np.int)

    # fields for empirical probabilities
    # counts
    nb = np.zeros((2, 4))
    nbxx = np.zeros((2, 4))

    # class boundaries for histograms
    # wclassbounds = np.linspace(0, 1, num=8)
    wlower = np.array([0,
                       0.1429,
                       0.2857,
                       0.4286,
                       0.5714,
                       0.7143,
                       0.8571])  # wclassbounds[0:7]
    wupper = np.array([0.1429,
                       0.2857,
                       0.4286,
                       0.5714,
                       0.7143,
                       0.8571,
                       1.0])  # wclassbounds[1:8]

    # evaluate mean rainfall intensity for wet boxes
    # these values should be determined during the aggregation phase!!!!!
    # mean volume threshold
    meanvol = np.percentile(vdn.precip[vdn.precip > 0.],
                            cascade_opt.percentile)  # np.mean(vdn.precip[vdn.precip>0.])
    cascade_opt.threshold = np.array([meanvol])

    # 2nd step: classify boxes at the upper level
    for i in range(0, n_out):
        if vdn.precip.values[i] > 0.:  # rain?
            if i == 0:  # only starting or isolated
                if vdn.precip.values[i+1] > 0.:
                    vbtype[i] = cascade.BoxTypes.starting
                else:
                    vbtype[i] = cascade.BoxTypes.isolated
            elif i == n_out-1:  # only ending or isolated
                if vdn.precip.values[i-1] > 0.:
                    vbtype[i] = cascade.BoxTypes.ending
                else:
                    vbtype[i] = cascade.BoxTypes.isolated
            else:  # neither at at the end nor at the beginning
                if vdn.precip.values[i-1] == 0. and vdn.precip.values[i+1] == 0.:
                    vbtype[i] = cascade.BoxTypes.isolated
                if vdn.precip.values[i-1] == 0. and vdn.precip.values[i+1] > 0.:
                    vbtype[i] = cascade.BoxTypes.starting
                if vdn.precip.values[i-1] > 0. and vdn.precip.values[i+1] > 0.:
                    vbtype[i] = cascade.BoxTypes.enclosed
                if vdn.precip.values[i-1] > 0. and vdn.precip.values[i+1] == 0.:
                    vbtype[i] = cascade.BoxTypes.ending
        else:
            vbtype[i] = cascade.BoxTypes.dry  # no rain

    # 3rd step: examine branching
    j = 0
    for i in range(0, n_in):
        if np.mod(i, 2) != 0:
            if vdn.precip.values[j] > 0:
                if vdn.precip.values[j] > meanvol:
                    belowabove = 1  # above mean
                else:
                    belowabove = 0  # below mean

                nb[belowabove, vbtype[j]-1] += 1

                if vec_data.precip.values[i-1] > 0 and vec_data.precip.values[i] == 0:
                    # P(1/0)
                    cascade_opt.p10[belowabove, vbtype[j]-1] += 1

                if vec_data.precip.values[i-1] == 0 and vec_data.precip.values[i] > 0:
                    # P(0/1)
                    cascade_opt.p01[belowabove, vbtype[j]-1] += 1

                if vec_data.precip.values[i-1] > 0 and vec_data.precip.values[i] > 0:
                    # P(x/x)
                    cascade_opt.pxx[belowabove, vbtype[j]-1] += 1

                    nbxx[belowabove, vbtype[j]-1] += 1

                    # weights
                    r1 = vec_data.precip.values[i-1]
                    r2 = vec_data.precip.values[i]
                    wxxval = r1 / (r1 + r2)

                    # Test
                    if abs(r1+r2-vdn.precip.values[j]) > 1.E-3:
                        print('i=' + str(i) + ', j=' + str(j) +
                              ', r1=' + str(r1) + ", r2=" + str(r2) +
                              ", Summe=" + str(vdn.precip.values[j]))
                        print(vec_data.index[i])
                        print(vdn.index[j])
                        print('error')
                        return cascade_opt, vdn

                    for k in range(0, 7):
                        if wxxval > wlower[k] and wxxval <= wupper[k]:
                            cascade_opt.wxx[k, belowabove, vbtype[j]-1] += 1
                            break
            j = j + 1

    # 4th step: transform counts to percentages
    cascade_opt.p01 = cascade_opt.p01 / nb
    cascade_opt.p10 = cascade_opt.p10 / nb
    cascade_opt.pxx = cascade_opt.pxx / nb

    with np.errstate(divide='ignore', invalid='ignore'):  # do not issue warnings here when dividing by zero, this is handled below
        for k in range(0, 7):
            cascade_opt.wxx[k, :, :] = cascade_opt.wxx[k, :, :] / nbxx[:, :]

    # In some cases, the time series are too short for deriving statistics.

    if (np.isnan(cascade_opt.p01).any() or
            np.isnan(cascade_opt.p10).any() or
            np.isnan(cascade_opt.pxx).any()):
        print("ERROR (branching probabilities):")
        print("Invalid statistics. Default values will be returned. "
              "Try to use longer time series or apply statistics "
              "derived for another station.")
        cascade_opt.fill_with_sample_data()

    # For some box types, the corresponding probabilities might yield nan.
    # If this happens, nan values will be replaced by 1/7 in order to provide
    # valid values for disaggregation.
    if np.isnan(cascade_opt.wxx).any():
        print("Warning (weighting probabilities):")
        print("The derived cascade statistics are not valid as some "
              "probabilities are undefined! ", end="")
        print("Try to use longer time series that might be more "
              "appropriate for deriving statistics. ", end="")
        print("As a workaround, default values according to equally "
              "distributed probabilities ", end="")
        print("will be applied...", end="")
        cascade_opt.wxx[np.isnan(cascade_opt.wxx)] = 1.0 / 7.0
        wxx = np.zeros((2, 4))
        for k in range(0, 7):
            wxx[:, :] += cascade_opt.wxx[k, :, :]
        if wxx.any() > 1.001 or wxx.any() < 0.999:
            print("failed! Using default values!")
            cascade_opt.fill_with_sample_data()
        else:
            print("OK!")

    return cascade_opt, vdn


def seasonal_subset(dataframe,
                    months='all'):
    '''Get the seasonal data.

    Parameters
    ----------
    dataframe : pd.DataFrame
    months: int, str
        Months to use for statistics, or 'all' for 1-12 (default='all')
    '''

    if isinstance(months, str) and months == 'all':
        months = np.arange(12) + 1

    for month_num, month in enumerate(months):
        df_cur = dataframe[dataframe.index.month == month]

        if month_num == 0:
            df = df_cur
        else:
            df = df.append(df_cur)

    return df.sort_index()


def build_casc(ObsData, hourly=True,level=9,
               months=None,
               avg_stats=True,
               percentile=50):
    '''Builds the cascade statistics of observed data for disaggregation

    Parameters
    -----------
    ObsData : pd.Series
        hourly=True -> hourly obs data
        else -> 5min data (disaggregation level=9 (default), 10, 11)
    
    months : numpy array of ints
        Months for each seasons to be used for statistics (array of
        numpy array, default=1-12, e.g., [np.arange(12) + 1])
    avg_stats : bool
        average statistics for all levels True/False (default=True)
    percentile : int, float
        percentile for splitting the dataset in small and high
        intensities (default=50)

    Returns
    -------
    list_seasonal_casc :
        list holding the results
    '''

    list_seasonal_casc = list()

    if months is None:
        months = [np.arange(12) + 1]

    # Parameter estimation for each season
    for cur_months in months:
        vdn = seasonal_subset(ObsData, cur_months)
        if len(ObsData.precip[np.isnan(ObsData.precip)]) > 0:
            ObsData.precip[np.isnan(ObsData.precip)] = 0

        casc_opt = melodist.cascade.CascadeStatistics()
        casc_opt.percentile = percentile
        list_casc_opt = list()
        
        count = 0
        
        if hourly:
            aggre_level = 5
        else:
            aggre_level =  level
        
        thresholds = np.zeros(aggre_level) #np.array([0., 0., 0., 0., 0.])
        
        for i in range(0, aggre_level):
            # aggregate the data
            casc_opt_i, vdn = aggregate_precipitation(vdn, hourly, \
                                percentile=percentile)
            thresholds[i] = casc_opt_i.threshold
            copy_of_casc_opt_i = copy.copy(casc_opt_i)
            list_casc_opt.append(copy_of_casc_opt_i)
            n_vdn = len(vdn)
            casc_opt_i * n_vdn  # level related weighting
            casc_opt + casc_opt_i  # add to total statistics
            count = count + n_vdn
        casc_opt * (1. / count)  # transfer weighted matrices to probabilities
        casc_opt.threshold = thresholds
        
        # statistics object
        if avg_stats:
            # in this case, the average statistics will be applied for all levels likewise
            stat_obj = casc_opt
        else:
            # for longer time series, separate statistics might be more appropriate
            # level dependent statistics will be assumed
            stat_obj = list_casc_opt

        list_seasonal_casc.append(stat_obj)

    return list_seasonal_casc

def sub_level_index(daily_index, level=9, fill_gaps=False):
    frequency = 42187500 * (2**(11-level))
    sublevel_index = pd.date_range(start=daily_index.min(), end=daily_index.max().replace(hour=23,minute=59,second=59), freq= '%sU'% frequency)
    # remove days that are not in the daily index:
    time_steps = 2**level
    if not fill_gaps and len(sublevel_index) > len(daily_index) * time_steps:
        df = pd.DataFrame(index=sublevel_index.date, data=dict(hour=sublevel_index.hour,\
                                                               minute=sublevel_index.minute,second=sublevel_index.second, keep=True))
        df.index = pd.to_datetime(df.index) #got a warning at firsr: "to_datetime"is depreciated
        dropdates = list(set(sublevel_index.date) - set(daily_index.date))
        df.loc[dropdates, 'keep'] = False
        df = df[df.keep]
        sublevel_index = pd.DatetimeIndex([pd.datetime(y, m, d, h,minute,second) for y, m, d, h,minute,second in zip(df.index.year,
                                                                                df.index.month,
                                                                                df.index.day,
                                                                                df.hour,df.minute,df.second)])

    return sublevel_index

def fmin_index(daily_index, fill_gaps=False):
    fminindex = pd.date_range(start=daily_index.min(), end=daily_index.max().replace(hour=23,minute=59,second=59), freq='300000L')
    # remove days that are not in the daily index:
    if not fill_gaps and len(fminindex) > len(daily_index) * 288:
        df = pd.DataFrame(index=fminindex.date, data=dict(hour=fminindex.hour,minute=fminindex.minute,second=fminindex.second, keep=True))
        df.index = pd.to_datetime(df.index) #got a warning at firsr: "to_datetime"is depreciated
        dropdates = list(set(fminindex.date) - set(daily_index.date))
        df.loc[dropdates, 'keep'] = False
        df = df[df.keep]
        fminindex = pd.DatetimeIndex([pd.datetime(y, m, d, h,minute,second) for y, m, d, h,minute,second in zip(df.index.year,
                                                                                df.index.month,
                                                                                df.index.day,
                                                                                df.hour,df.minute,df.second)])

    return fminindex

