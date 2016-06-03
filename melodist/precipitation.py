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
import copy
import melodist
import melodist.util
import numpy as np
import pandas as pd

def disagg_prec(dailyData, method='equal', cascade_options=None, hourly_data_obs=None, zerodiv="uniform", shift=0):
    """The disaggregation function for precipitation
    
    Parameters
    ----
    dailyData :       daily data     
    method :          method to dissaggregate 
    cascade_options*: cascade object including statistical parameters for the cascade model
    hourly_data_obs*: observed hourly data of master station
    zerodiv*:         method to deal with zero division by key "uniform" --> uniform distribution
    shift*:           shifts the precip data by shift (int) steps (eg +7 for 7:00 to 6:00)    
    """

    if method not in ('equal', 'cascade', 'masterstation'):
        raise ValueError('Invalid option')

    if method == 'equal':
        precip_disagg = melodist.distribute_equally(dailyData.precip, divide=True)
    elif method == 'masterstation':
        precip_disagg = precip_master_station(dailyData, hourly_data_obs, zerodiv)       
    elif method == 'cascade':
        assert cascade_options is not None
        precip_disagg = disagg_prec_cascade(dailyData, cascade_options, shift=shift)

    return precip_disagg

def disagg_prec_cascade(precip_daily, cascade_options,shift=0, test=False):
    """Precipitation disaggregation using the cascade model proposed by
    Olsson (1998)
    
    Parameters
    ----
    precip_daily :    daily data        
    cascade_options*: cascade object including statistical parameters for the cascade model
    shift:            shifts the precip data by shift (int) steps (eg +7 for 7:00 to 6:00)
    test:             test mode, returns time series of each cascade level
    """

    if len(precip_daily)  < 2:
        raise ValueError('Input data must have at least two elements.')

    # set missing values to zero:
    precip_daily = precip_daily.copy()
    missing_days = precip_daily.index[precip_daily.isnull()]
    precip_daily[missing_days] = 0
    
    si = 5 # index of first level    
    
    # statistics for branching into two bins
    wxxcum = np.zeros((7,2,4))
    
    if(isinstance(cascade_options, melodist.cascade.CascadeStatistics)):
        # this is the standard case considering one data set for all levels
        # get cumulative probabilities for branching
        overwrite_stats = False
        for k in range(0,7):
            wxxcum[k,:,:] = cascade_options.wxx[k,:,:]
            if k > 0:
                wxxcum[k,:,:] = wxxcum[k-1,:,:] + wxxcum[k,:,:];        
    elif(isinstance(cascade_options, list)):     
        if(len(cascade_options) == 5):
            overwrite_stats = True
            list_casc = cascade_options
        else:
            raise ValueError('Cascade statistics list must have 5 elements!')
    else:
        raise TypeError('cascade_options has invalid type')
    
    # arrays for each level
    n = len(precip_daily)
    vdn1 = np.zeros(n*2)
    vdn2 = np.zeros(n*4)
    vdn3 = np.zeros(n*8)
    vdn4 = np.zeros(n*16)
    vdn5 = np.zeros(n*32)

    # class boundaries for histograms
    wclassbounds = np.array([0.0,0.1429,0.2857,0.4286,0.5714,0.7143,0.8571,1.0]) #np.linspace(0,1,num=8)

    # disaggregation for each level
    for l in range(1,6):
        if l == 1:
            vdn_in  = precip_daily
            vdn_out = vdn1;
        elif l == 2:
            vdn_in  = vdn_out;
            vdn_out = vdn2;
        elif l == 3:
            vdn_in  = vdn_out;
            vdn_out = vdn3;
        elif l == 4:
            vdn_in  = vdn_out;
            vdn_out = vdn4;
        elif l == 5:
            vdn_in  = vdn_out;
            vdn_out = vdn5;
        
        si -= 1
        if(overwrite_stats):
            cascade_options = list_casc[si]
            for k in range(0,7):
                wxxcum[k,:,:] = cascade_options.wxx[k,:,:]
                if k > 0:
                    wxxcum[k,:,:] = wxxcum[k-1,:,:] + wxxcum[k,:,:]
            meanvol = cascade_options.threshold[0]
        else:
            meanvol = cascade_options.threshold[si]
        
        # evaluate time step
        dt = (precip_daily.index[1] - precip_daily.index[0]).total_seconds() / 3600 # hours
        print("disaggregating " + str(dt/2**(l-1)) + " hours to " + str(dt/2**l) + " hours...")

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
                if i==0: # only starting or isolated
                    if vdn_in[i+1] > 0:
                        vbtype = cascade.BoxTypes.starting
                    else:
                        vbtype = cascade.BoxTypes.isolated
                elif i==len(vdn_in)-1: # only ending or isolated
                    if vdn_in[i-1] > 0:
                        vbtype = cascade.BoxTypes.ending
                    else:
                        vbtype = cascade.BoxTypes.isolated
                else: # neither at at the end nor at the beginning
                    if vdn_in[i-1] == 0 and vdn_in[i+1] == 0:
                        vbtype = cascade.BoxTypes.isolated

                    if vdn_in[i-1] == 0 and vdn_in[i+1] > 0:
                        vbtype = cascade.BoxTypes.starting

                    if vdn_in[i-1] > 0 and vdn_in[i+1] > 0:
                        vbtype = cascade.BoxTypes.enclosed;

                    if vdn_in[i-1] > 0 and vdn_in[i+1] == 0:
                        vbtype = cascade.BoxTypes.ending;

                # above or below mean?
                if vdn_in[i] > meanvol:
                    belowabove = 1 # above mean
                else:
                    belowabove = 0 # below mean

                #
                p = np.zeros((3,1))
                p[0] = cascade_options.p01[belowabove,vbtype-1] # index changed!
                p[1] = cascade_options.p10[belowabove,vbtype-1]
                p[2] = cascade_options.pxx[belowabove,vbtype-1]

                # draw a random number to determine the braching type
                rndp = np.random.random()

                if rndp <= p[0]:
                    # first box 0, second box: 1 P(0/1)
                    vdn_out[j] = 0.0
                    j = j + 1
                    vdn_out[j] = vdn_in[i]
                    j = j + 1
                elif rndp > p[0] and rndp <= p[0] + p[1]: # "p[0] +" added!!!!!!!!!!!!!!!!!!
                    # first box 1, second box: 0 P(1/0)
                    vdn_out[j] = vdn_in[i]
                    j = j + 1
                    vdn_out[j] = 0.0
                    j = j +1
                else:
                    # both boxes wet
                    # we need a new random number
                    rndw = np.random.random()
                    
                    # guess w1:
                    for k in range(0,7):
                        if rndw <= wxxcum[k,belowabove,vbtype-1]:
                            w1 = wclassbounds[k+1] - 1./14. # class center
                            break

                    vdn_out[j] = w1 * vdn_in[i]
                    j = j + 1;
                    vdn_out[j] = (1. - w1) * vdn_in[i]
                    j = j + 1

                    # check results (in the previous version this error has never been observed)
                    if w1 < 0 or w1 > 1:
                        print('error')
                        return
            else:
                # add two dry boxes
                vdn_out[j] = 0.0
                j = j + 1;
                vdn_out[j] = 0.0
                j = j + 1;


    # uniformly disaggregate 0.75 h values to 0.25 h values
    print('disaggregating 0.75 hours to 0.25 hours... (uniform)')
    vdn_025 = np.zeros(len(vdn_out)*3);
    # vtn_025 = zeros(length(vdn_out)*3,1);
    # dt = (vtn_out(2) - vtn_out(1))/3;
    j = 0;
    for i in range (0,len(vdn_out)):
        for m in range(0,3):
            vdn_025[j+m] = vdn_out[i] / 3.;
            #vtn_025(j+m) = vtn_out(i) - 2 * dt + m *dt;
        j = j + 3;
    
    # aggregate to hourly time steps
    print('aggregating 0.25 hours to 1.0 hours...')
    vdn_025cs = np.cumsum(vdn_025)
    vdn = np.zeros(int(len(vdn_025)/4))
    for i in range(0,len(vdn)+1):
        #for first hour take 4th item
        if i==0:  
            vdn[i] = vdn_025cs[3]
        #pass 1
        elif i==1:
            pass
        
        else: 
            #>1 (starting with 2-1 = 1 item)
            vdn[i-1] = vdn_025cs[(i*4)-1] - vdn_025cs[(i*4)-5]

    precip_hourly = pd.Series(index=melodist.util.hourly_index(precip_daily.index), data=vdn)

    # set missing days to nan again:
    for date in missing_days:
        precip_hourly[precip_hourly.index.date == date.date()] = np.nan
    
    #shifts the data by shift steps (fills with nan/cuts edge data )
    if shift != 0:
        precip_hourly = precip_hourly.shift(shift)

    # replace column
    print("done")
    
    # return time series
    if test:
        return vdn1, vdn2, vdn3, vdn4, vdn5, vdn_025, precip_hourly
    else:
        return precip_hourly

###############################################################################
#disaggregate by master station
def precip_master_station(precip_daily, master_precip_hourly, zerodiv):
    """Disaggregate precipitation based on the patterns of a master station

    Parameters
    ----
    precip_daily :   daily data 
    master_precip_hourly :  observed hourly data of the master station
    zerodiv :        method to deal with zero division by key "uniform" --> uniform distribution
    """

    precip_hourly = pd.Series(index=melodist.util.hourly_index(precip_daily.index))
    
    # set some parameters for cosine function
    for index_d, precip in precip_daily.iteritems():
        
        #get hourly data of the day
        index = index_d.date().isoformat()
        precip_h = master_precip_hourly[index]
        
        #calc rel values and multiply by daily sums    
        #check for zero division    
        if precip_h.sum() != 0 and precip_h.sum() != np.isnan(precip_h.sum()):            
            precip_h_rel = (precip_h / precip_h.sum()) * precip
            
        else:
            #uniform option will preserve daily data by uniform distr
            if zerodiv == 'uniform':
                precip_h_rel = (1/24) * precip
    
            else:
                precip_h_rel = 0
                       
        #write the disaggregated day to data
        precip_hourly[index] = precip_h_rel

    return precip_hourly


def aggregate_precipitation(vec_data):
    """Aggregates highly resolved precipitation data and creates statistics
    
        Args:
        vec_data: hourly values

    Returns:
        cascade object representing statistics of the cascade model
    """
    cascade_opt = cascade.CascadeStatistics()

    # get time step
    #dt = vec_data.index[0].offset.nanos / 1E9 / 3600 # hours
    dt = (vec_data.index[1]-vec_data.index[0]).total_seconds() / 3600 # hours

    dt_out = dt * 2

    print('aggregating ' + str(dt) + ' hours to ' + str(dt_out) + ' hours...')

    # length of input time series
    n_in = len(vec_data)
    n_out = np.floor(n_in/2);

    # alternative:
    # 1st step: new time series
    vec_time = vec_data.index
    vdn0=[]
    vtn0=[]
    j = 0
    for i in range(0,n_in):
        if(np.mod(i,2)!=0):
            vdn0.append(vec_data.precip[i-1] + vec_data.precip[i])
            vtn0.append(vec_time[i])
            j=j+1

    vdn = pd.DataFrame(index=vtn0, data={'precip':vdn0})

    # length of new time series
    n_out = len(vdn)

    #print('n_in='+str(n_in)+', n_out=' + str(n_out))

    # series of box types:
    vbtype = np.zeros((n_out,), dtype=np.int)

    # fields for empirical probabilities
    # counts
    nb  = np.zeros((2,4))
    nbxx  = np.zeros((2,4))

    # class boundaries for histograms
    # wclassbounds = np.linspace(0,1,num=8)
    wlower = np.array([0,0.1429,0.2857,0.4286,0.5714,0.7143,0.8571]) # wclassbounds[0:7]
    wupper = np.array([0.1429,0.2857,0.4286,0.5714,0.7143,0.8571,1.0]) # wclassbounds[1:8]

    # evaluate mean rainfall intensity for wet boxes
    # these values should be determined during the aggregation phase!!!!!
    # mean volume threshold
    meanvol = np.percentile(vdn.precip[vdn.precip>0.], cascade_opt.percentile) # np.mean(vdn.precip[vdn.precip>0.])
    cascade_opt.threshold = np.array([meanvol])

    # 2nd step: classify boxes at the upper level
    for i in range(0, n_out):
        if vdn.precip[i] > 0.: # rain?
            if i == 0: # only starting or isolated
                if vdn.precip[i+1] > 0.:
                    vbtype[i] = cascade.BoxTypes.starting
                else:
                    vbtype[i] = cascade.BoxTypes.isolated
            elif i == n_out-1: # only ending or isolated
                if vdn.precip[i-1] > 0.:
                    vbtype[i] = cascade.BoxTypes.ending
                else:
                    vbtype[i] = cascade.BoxTypes.isolated
            else: # neither at at the end nor at the beginning
                if vdn.precip[i-1] == 0. and vdn.precip[i+1] == 0.:
                    vbtype[i] = cascade.BoxTypes.isolated
                if vdn.precip[i-1] == 0. and vdn.precip[i+1] > 0.:
                    vbtype[i] = cascade.BoxTypes.starting
                if vdn.precip[i-1] > 0. and vdn.precip[i+1] > 0.:
                    vbtype[i] = cascade.BoxTypes.enclosed
                if vdn.precip[i-1] > 0. and vdn.precip[i+1] == 0.:
                    vbtype[i] = cascade.BoxTypes.ending
        else:
            vbtype[i] = cascade.BoxTypes.dry   # no rain

    # 3rd step: examine branching
    j=0
    for i in range(0, n_in):
        if np.mod(i,2) != 0:
            if vdn.precip[j] > 0:
                if vdn.precip[j] > meanvol:
                    belowabove = 1 # above mean
                else:
                    belowabove = 0 # below mean

                nb[belowabove,vbtype[j]-1] += 1

                if vec_data.precip[i-1] > 0 and vec_data.precip[i] == 0:
                    # P(1/0)
                    cascade_opt.p10[belowabove,vbtype[j]-1] += 1

                if vec_data.precip[i-1] == 0 and vec_data.precip[i] > 0:
                    # P(0/1)
                    cascade_opt.p01[belowabove,vbtype[j]-1] += 1

                if vec_data.precip[i-1] > 0 and vec_data.precip[i] > 0:
                    # P(x/x)
                    cascade_opt.pxx[belowabove,vbtype[j]-1] += 1

                    nbxx[belowabove,vbtype[j]-1] += 1

                    # weights
                    r1 = vec_data.precip[i-1]
                    r2 = vec_data.precip[i]
                    wxxval = r1 / (r1 + r2)
                    #print(wxxval)

                    # Test
                    if abs(r1+r2-vdn.precip[j]) > 1.E-3:
                        print('i='+str(i)+', j='+str(j)+', r1='+str(r1)+", r2="+str(r2)+", Summe="+str(vdn.precip[j]))
                        print(vec_data.index[i])
                        print(vdn.index[j])
                        print('error')
                        return cascade_opt, vdn

                    for k in range(0,7):
                        if wxxval > wlower[k] and wxxval <= wupper[k]:
                            cascade_opt.wxx[k,belowabove,vbtype[j]-1] += 1
                            break
            j = j + 1

    # 4th step: transform counts to percentages
    cascade_opt.p01 = cascade_opt.p01 / nb;
    cascade_opt.p10 = cascade_opt.p10 / nb;
    cascade_opt.pxx = cascade_opt.pxx / nb;

    for k in range(0,7):
        cascade_opt.wxx[k,:,:] = cascade_opt.wxx[k,:,:] / nbxx[:,:]

    # In some cases, the time series are too short for deriving statistics.

    if np.isnan(cascade_opt.p01).any() or np.isnan(cascade_opt.p10).any() or np.isnan(cascade_opt.pxx).any():
        print("ERROR (branching probabilities):")
        print("Invalid statistics. Default values will be returned. Try to use longer time series or apply statistics derived for another station.")
        cascade_opt.fill_with_sample_data()

    # For some box types, the corresponding probabilities might yield nan.
    # If this happens, nan values will be replaced by 1/7 in order to provide
    # valid values for disaggregation.
    if np.isnan(cascade_opt.wxx).any():
        print("Warning (weighting probabilities):")
        print("The derived cascade statistics are not valid as some probabilities are undefined! ", end="")
        print("Try to use longer time series that might be more appropriate for deriving statistics. ", end="")
        print("As a workaround, default values according to equally distributed probabilities ", end="")
        print("will be applied...", end="")
        cascade_opt.wxx[np.isnan(cascade_opt.wxx)] = 1.0 / 7.0
        wxx = np.zeros((2,4))
        for k in range(0,7):
            wxx[:,:] += cascade_opt.wxx[k,:,:]
        if wxx.any() > 1.001 or wxx.any() < 0.999:
            print("failed! Using default values!")
            cascade_opt.fill_with_sample_data()
        else:
            print("OK!")

    return cascade_opt, vdn


def seasonal_subset(dataframe, months='all'):
    '''get the data seasonal

    Parameters
    ----
    dataframe:
    months:    Months to use for statistics, or 'all' for 1-12 (default='all')

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


def build_casc(hourlyDataObs, months=None, avg_stats=True, percentile=50):
    '''builds the cascade satistics of the observed data for fallowing disaggregation

    Parameters
    ----
    hourlyDataObs : hourly obserevd data
    months :        Months for each seasons to be used for statistics (array of numpy array, default=1-12, e.g., [np.arange(12) + 1])
    avg_stats :     average statistics for all levels True/False (default=True)
    percentile :    percentil for splitting the dataset in small and high intensities (default=50)

    Returns
    ----
    list_seasonal_casc : list holding the results
    '''

    list_seasonal_casc = list()

    if months is None:
        months = [np.arange(12) + 1]

    # Parameter estimation for each season
    for cur_months in months:
        vdn = seasonal_subset(hourlyDataObs,cur_months)
        if(len(hourlyDataObs.precip[np.isnan(hourlyDataObs.precip)]) > 0):
            hourlyDataObs.precip[np.isnan(hourlyDataObs.precip)]=0

        casc_opt = melodist.cascade.CascadeStatistics()
        casc_opt.percentile = percentile
        list_casc_opt = list()
        thresholds = np.array([0., 0., 0., 0., 0.])
        count = 0
        for i in range(0,5):
            #aggregate the data
            casc_opt_i,vdn = aggregate_precipitation(vdn)
            thresholds[i] = casc_opt_i.threshold
            copy_of_casc_opt_i = copy.copy(casc_opt_i)
            list_casc_opt.append(copy_of_casc_opt_i)
            n_vdn = len(vdn)
            casc_opt_i * n_vdn # level related weighting
            casc_opt + casc_opt_i # add to total statistics
            count = count + n_vdn
        casc_opt * (1. / count) # transfer weighted matrices to probabilities
        casc_opt.threshold = thresholds

        # statistics object
        if(avg_stats):
            # in this case, the average statistics will be applied for all levels likewise
            stat_obj = casc_opt
        else:
            # for longer time series, separate statistics might be more appropriate
            # level dependent statistics will be assumed
            stat_obj = list_casc_opt

        list_seasonal_casc.append(stat_obj)

    return list_seasonal_casc
