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

import melodist
import pandas as pd
import numpy as np
import collections

###############################################################################
#reads smet data

def read_smet(filename, mode):
    """Reads smet data and returns the data in required dataformat (pd df)
    
    See https://models.slf.ch/docserver/meteoio/SMET_specifications.pdf
    for further details on the specifications of this file format.

    Parameters
    ----
    filename : SMET file to read
    mode :     "d" for daily and "h" for hourly input

    Returns
    ----
    [header, data]
    header:    header as dict
    data :     data as pd df
    """

    # dictionary
    #based on smet spec V.1.1 and self defined
    #daily data
    dict_d=   {'TA':'tmean',
               'TMAX':'tmax',   #no spec
               'TMIN':'tmin',   #no spec
               'PSUM':'precip',
               'ISWR':'glob',     #no spec
               'RH':'hum',
               'VW':'wind'
                }

    #hourly data
    dict_h=   {'TA':'temp',
               'PSUM':'precip',
               'ISWR':'glob',     #no spec
               'RH':'hum',
               'VW':'wind'
                }

    with open(filename) as f:
        in_header = False
        data_start = None
        header = collections.OrderedDict()

        for line_num, line in enumerate(f):

            if line.strip() == '[HEADER]':
                in_header = True
                continue
            elif line.strip() == '[DATA]':
                data_start = line_num + 1
                break

            if in_header:
                line_split = line.split('=')
                k = line_split[0].strip()
                v = line_split[1].strip()
                header[k] = v


    #get column names
    columns = header['fields'].split()
    multiplier = [float(x) for x in header['units_multiplier'].split()][1:]
    

    data = pd.read_table(
        filename,
        sep=r'\s+',
        na_values=[-999],
        skiprows=data_start,
        names=columns,
        index_col='timestamp',
        parse_dates=True,
        )

    data = data*multiplier

    del data.index.name

    #rename columns
    if mode == "d":
        data = data.rename(columns=dict_d)
    if mode == "h":
        data = data.rename(columns=dict_h)

    return header, data
    
################################################################################
#read dwd data

def read_dwd(filename, metadata, mode="d", skip_last=True):
    """Reads dwd (German Weather Service) data and returns the data in required
    dataformat (pd df)

    Parameters
    ----
    filename : DWD file to read (full path) / list of hourly files (RR+TU+FF) 
    metadata : corresponding DWD metadata file to read
    mode :    "d" for daily and "h" for hourly input
    skip_last : boolen, skips last line due to file format

    Returns
    ----
    [header, data]
    header:    header as dict
    data :     data as pd df
    """
    
    #--------------------------------------------------------------------------
    def read_single_dwd(filename, metadata, mode, skip_last):

        #Param names {'DWD':'dissag_def'}
        dict_d=   {'LUFTTEMPERATUR':'tmean',
                   'LUFTTEMPERATUR_MINIMUM':'tmin',   #no spec
                   'LUFTTEMPERATUR_MAXIMUM':'tmax',   #no spec
                   'NIEDERSCHLAGSHOEHE':'precip',
                   'XXXXXXX':'glob',     #no spec
                   'REL_FEUCHTE':'hum',
                   'WINDGESCHWINDIGKEIT':'wind',
                   'SONNENSCHEINDAUER':'sun_h'
                    }
        
        #---read meta------------------
        meta = pd.read_csv(
            metadata,
            sep=';'
            )
        
        #remove whitespace from header columns
        meta.rename(columns=lambda x: x.strip(), inplace=True)
        
        header = {"Stations_id" : meta.Stations_id[meta.last_valid_index()],
                  "Stationsname": meta.Stationsname[meta.last_valid_index()],
                  #workaround for colnames with . (Geogr.Breite)
                  "Breite" : meta.iloc[meta.last_valid_index(),2], #DezDeg
                  "Laenge" : meta.iloc[meta.last_valid_index(),3]  #DezDeg
                  }
        
        #---read data------------------
        if skip_last!=None:
            num_lines = sum(1 for line in open(filename))
            skip_last = [num_lines-1]
        
        #hourly data must be parsed by custom definition
        if mode == "d":
            data = pd.read_csv(
                filename,
                sep=';',
                na_values= '-999',
                index_col=' MESS_DATUM',
                parse_dates=True,
                skiprows=skip_last
                )
        
        #hourly data must be parsed by custom definition
        if mode == "h":
    
            def date_parser(date_time):
                hour = date_time[8:10]
                day = date_time[6:8]
                month= date_time[4:6]
                year = date_time[0:4]
                minute = '00'
                sec = '00'
                return  pd.Timestamp('%s-%s-%s %s:%s:%s' % (year, month, day, hour, minute, sec))     
            
            data = pd.read_csv(
                filename,
                sep=';',
                na_values= '-999',
                index_col=' MESS_DATUM',
                date_parser=date_parser,
                skiprows=skip_last
                )
        
        #remove whitespace from header columns
        data.rename(columns=lambda x: x.strip(), inplace=True)
        
        #rename to dissag definition
        data = data.rename(columns=dict_d)
        #get colums which are not defined 
        drop = [col for col in data.columns if col not in dict_d.values()]
        #delete columns
        data = data.drop(drop, axis=1)
        
        
#        #convert temperatures to Kelvin (+273.15)
        if 'tmin' in data.columns:
            data["tmin"] = data["tmin"] + 273.15
        if 'tmax' in data.columns:
            data["tmax"] = data["tmax"] + 273.15
        if 'tmean' in data.columns:
            data["tmean"] = data["tmean"] + 273.15
        if 'temp' in data.columns:
            data["temp"] = data["temp"] + 273.15
        
        return header, data
    
    #-------------------------------------------------------------------------- 
    
    if type(filename)==list:
        i = 1
        for file in filename:
            header, data_h = read_single_dwd(file, metadata, mode, skip_last)
            
            if i == 1:
                data = data_h           
            else:          
                data = data.join(data_h, how='outer')
            i+=1
                
    else:
        header, data = read_single_dwd(filename, metadata, mode, skip_last)
        
    return header, data

################################################################################
#reads csv data

def read_csv(filename, mode):
    """Reads smet data and returns the data in required dataformat (pd df)

    Parameters
    ----
    filename : SMET file to read
    mode :     "d" for daily and "h" for hourly input

    Returns
    ----
    data :     data as pd df
    """

    if mode == "d":

        daily_data_raw = pd.read_csv(
            filename,
            index_col='datum',
            parse_dates=True,
            na_values=[9999, 32700],
            )

        del daily_data_raw.index.name

        daily_data = pd.DataFrame(
            index=daily_data_raw.index,
            columns=melodist.columns_daily,
            data=dict(
                tmin=daily_data_raw.tmin / 10,
                tmax=daily_data_raw.tmax / 10,
                tmean=daily_data_raw.t / 10,
                precip=daily_data_raw.nied.replace(-1, 0) / 10,
                glob=daily_data_raw.strahl * 10000 / 3600 / 24., # radiation is in J/cm2, convert to mean hourly radiation in W/m2
                hum=daily_data_raw.rel,
                wind=daily_data_raw.vv.replace(999, np.nan) / 10,
                ),
            )

        daily_data[['tmin', 'tmax', 'tmean']] += 273.15 # convert temperatures to K

        return daily_data


    #hourly mode
    if mode == "h":
        hourly_data_obs = pd.read_csv(
            filename,
            parse_dates=True,
            index_col='date',
        )

        del hourly_data_obs.index.name

        return hourly_data_obs


################################################################################
#writes smet-file
#time can hard-coded overwrite the times to day values "T00:00"

def write_smet(filename, data, metadata, nodata_value=-999, mode='h', check_nan=True):
    """writes smet files
    
    Parameters
    ----
    filename :    filename/loction of output
    data :        data to write as pandas df
    metadata:     header to write input as dict
    nodata_value: Nodata Value to write/use
    mode:         defines if to write daily ("d") or continuos data (default 'h') 
    check_nan:    will check if only nans in data and if true will not write this colums (default True)
    """
    
    # dictionary
    #based on smet spec V.1.1 and selfdefined
    #daily data
    dict_d=   {'tmean':'TA',
               'tmin':'TMAX',   #no spec
               'tmax':'TMIN',   #no spec
               'precip':'PSUM',
               'glob':'ISWR',     #no spec
               'hum':'RH',
               'wind':'VW'
                }

    #hourly data
    dict_h=   {'temp':'TA',
               'precip':'PSUM',
               'glob':'ISWR',     #no spec
               'hum':'RH',
               'wind':'VW'
                }
                
    #rename columns
    if mode == "d":
        data = data.rename(columns=dict_d)
    if mode == "h":
        data = data.rename(columns=dict_h)

    if check_nan:     
        #get all colums with data
        datas_in = data.sum().dropna().to_frame().T
        #get colums with no datas
        drop = [data_nan for data_nan in data.columns if data_nan not in datas_in]    
        #delete columns
        data = data.drop(drop, axis=1)
    
    with open(filename, 'w') as f:

        #preparing data
        #converte date_times to SMET timestamps
        if mode == "d":
            t = '%Y-%m-%dT00:00'
        if mode == "h":
            t = '%Y-%m-%dT%H:%M'

        data['timestamp'] = [d.strftime(t) for d in data.index]
        
        cols = data.columns.tolist()
        cols = cols[-1:] + cols[:-1]
        data = data[cols]


        #metadatas update
        metadata['fields'] = ' '.join(data.columns)
        metadata["units_multiplier"] = len(metadata['fields'].split())*"1 "

        #writing data
        #metadata
        f.write('SMET 1.1 ASCII\n')
        f.write('[HEADER]\n')

        for k, v in metadata.items():
            f.write('{} = {}\n'.format(k, v))

        #data
        f.write('[DATA]\n')

        data_str = data.fillna(nodata_value).to_string(
            header=False,
            index=False,
            float_format=lambda x: '{:.2f}'.format(x),
        )

        f.write(data_str)

###############################################################################
