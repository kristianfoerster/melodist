[![Build Status](https://github.com/kristianfoerster/melodist/actions/workflows/ci.yml/badge.svg)](https://github.com/kristianfoerster/melodist/actions/workflows/ci.yml)

# MELODIST - An open-source MEteoroLOgical observation time series DISaggregation Tool

*Welcome to MELODIST - MEteoroLOgical observation time series DISaggregation Tool*

MELODIST is an open-source toolbox written in Python for disaggregating daily meteorological time series to hourly time steps. It is licensed under GPLv3 (see license file). The software framework consists of disaggregation functions for each variable including temperature, humidity, precipitation, shortwave radiation, and wind speed. These functions can simply be called from a station object, which includes all relevant information about site characteristics. The data management of time series is handled using data frame objects as defined in the pandas package. In this way, input and output data can be easily prepared and processed. For instance, the pandas package is data i/o capable and includes functions to plot time series using the matplotlib library.

An [example file](examples/examples.ipynb) is provided along the package itself as a Jupyter notebook. This example demonstrates the usage of MELODIST for all variables. First, a station object is created providing some basic details on the site’s characteristics (e.g., latitude and longitude are relevant for radiation disaggregation). Once the station object is defined and filled with data, each disaggregation step is done through calling the designated function specified for each variable. Each of these functions requires a `method` argument and if needed additional parameters to work properly. Some of these methods (see below) require additional statistical evaluations of hourly time series prior to the disaggregation procedure. This information is stored in a station statistics object that is associated to the station object (see example file for further details).

## Station object
In the framework of MELODIST a station object includes all relevant information including metadata and time series. A station is generated using the constructor method:

```python
s = melodist.Station(lon=longitude, lat=latitude, timezone=timezone)
```

Data is simply added by assignment (e.g., the data frame `data_obs_daily`):

```python
s.data_daily = data_obs_daily
```

A station statistics object can be generated in a similar manner. As station statistics are derived through analysing hourly observations for calibration, a reference to the data frame including hourly observations is given:

```python
s.statistics = melodist.StationStatistics(data_obs_hourly)
```

Statistics can be derived for each variable by calling the respective function of the statistics object `s.statistics`: `calc_wind_stats()`, `calc_humidity_stats()`, `calc_temperature_stats()`, `calc_precipitation_stats()`, and `calc_radiation_stats`.

## Naming convention for dataframe columns

MELODIST expects exact naming conventions for the time series provided in pandas dataframes. Please find the specification of column names below:

* `temp`: Temperature [K]
* `precip`: Precipitation [mm/time step]
* `glob`: Global (shortwave) radiation [W/m<sup>2</sup>]
* `hum`: Relative humidity [%]
* `wind`: Wind speed [m/s]
* `ssd`: Sunshine duration [min]

For daily data, additional columns need to be specified (if applicable):

* `temp`: Average temperature [K]
* `tmin`: Minimum temperature [K]
* `tmax`: Maximum temperature [K]
* `hum`: Average humidity [%]
* `hum_min`: Minimum humidity [%]
* `hum_max`: Maximum humidity [%]
* `ssd`: Sunshine duration [h]

Please note that the dataframe's index must contain datetime values.

## Disaggregation methods

The `Station` class provides functions to perform the disaggregation procedures for each variable: `disaggregate_temperature()`, `disaggregate_humidity()`, `disaggregate_wind()`, `disaggregate_radiation()`, and `disaggregate_precipitation()`. Moreover, an interpolation approach is also available using the `interpolate()` function.

Hint: It is worth noting that each of the implemented disaggregation methods is directly accessible, e.g., `melodist.precipitation.disagg_prec()`. In this case all relevant parameters (e.g., those derived through calibrations) need to be provided in the function call. This method-specific call of functions is not necessary if a station and the corresponding station statistics object is defined. Thus, it is recommended to define objects and to perform the disaggregation procedures using the object’s methods. Also, the names and signatures of these functions are likely subject to changes in future versions of MELODIST.

Please find below a list of available disaggregation methods for each variable which can be specified in the respective disaggregation methods of a `Station` object:

### Temperature
* `method='sine_min_max'` (T1): standard sine redistribution; preserves T<sub>min</sub> and T<sub>max</sub> but not T<sub>mean</sub>.
* `method='sine_mean'`: sine redistribution; preserves T<sub>mean</sub> and the diurnal temperature range (T<sub>max</sub> – T<sub>min</sub>) but not T<sub>min</sub> and T<sub>max</sub>.
* `method='mean_course_min_max'`: redistribute following a prescribed temperature course calculated from hourly observations; preserves T<sub>min</sub> and T<sub>max</sub>.
* `method='mean_course_mean'`: redistribute following a prescribed temperature course calculated from hourly observations; preserves T<sub>mean</sub> and the diurnal temperature range.
* Possible options for `min_max_time` are:
  * `'fix'` (T1a): The diurnal course of temperature is fixed without any seasonal variations.
  * `'sun_loc'` (T1b): The diurnal course of temperature is modelled based on sunrise, noon and sunset calculations.
  * `'sun_loc_shift'` (T1c): This option activates empirical corrections of the ideal course modelled by `sun_loc` (requires calling `calc_temperature_stats()` prior to the disaggregation).
* An optional parameter `mod_nighttime` (T1d, bool, default: `False`) allows one to apply a linear interpolation of night time values, which proves preferable during polar nights.

### Humidity
* `method='equal'`: duplicate mean daily humidity for the 24 hours of the day.
* `method='minimal'` (H1): The dew point temperature is set to the minimum temperature on that day.
* `method='dewpoint_regression'` (H2): Based on hourly observations, a regression approach is applied to calculate daily dew point temperature. Regression parameters must be specified (which is automatically done if `calc_humidity_stats()` is called prior to disaggregation).
* `method='linear_dewpoint_variation'` (H3): This method extends H2 through linearly varying dew point temperature between consecutive days. The parameter `kr` needs to be specified (`kr=6` if monthly radiation exceeds 100 W/m<sup>2</sup> else `kr=12`).
* `method='min_max'` (H4): this method requires minimum and maximum relative humidity for each day.
* `method='month_hour_precip_mean'`: calculate hourly humidity from categorical [month, hour, precip(y/n)] mean values derived from observations.

### Wind speed
* `method='equal'` (W1): If this method is chosen, the daily average wind speed is assumed to be valid for each hour on that day.
* `method='cosine'` (W2): The cosine function option simulates a diurnal course of wind speed and requires calibration (`calc_wind_stats()`).
* `method='random'` (W3): This option is a stochastic method that draws random numbers to disaggregate wind speed taking into account the daily average (no parameter estimation required).

### Radiation
* `method='pot_rad'` (R1): This method allows one to disaggregate daily averages of shortwave radiation using hourly values of potential (clear-sky) radiation calculated for the location of the station.
* `method='pot_rad_via_ssd'` (R2): If daily sunshine recordings are available, the Angstrom model is applied to transform sunshine duration to shortwave radiation.
* `method='pot_rad_via_bc'` (R3): In this case, the Bristow-Campbell model is applied which relates minimum and maximum temperature to shortwave radiation.
* `method='mean_course'`: hourly radiation follows an observed average course (calculated for each month) while preserving the daily mean.

### Precipitation
* `method='equal'` (P1): In order to derive hourly from daily values, the daily total is simply divided by 24 resulting in an equal distribution.
* `method='cascade'` (P2): The cascade model is more complex and requires a parameter estimation method (`calc_precipitation_stats()`). Statistics can be calculated using different options (parameters). Using the keyword `months`, the seasons for which the statistics will be calculated independently can be specified (see example file). The keyword `percentile` allows one to adjust the threshold to separate precipitation intensities into two classes (low and high) for building the parameters. The default value is 50% (median). An additional optional argument `avg_stats` is used to decide whether statistics of all cascade levels will be averaged (default is `True`). All options previously listed are optional and can be changed to tune the disaggregation results. A new feature also allows for working with 5 minutes precipitation data. Both functions `disagg_prec_cascade()`and `aggregate_precipitation()` offer switching to 5min data by setting ``hourly=False`` in the function calls. Moreover, the former function allows for different level configurations (``levels`` can be set to 9 (standard), 10, or 11, depending on the number of branching levels to be considered in the cascade model). The usage of the cascade model for sub-hourly precipitation is demonstrated in a Jupyter notebook (see `/examples/precip5min_example.ipynb`).
* `method='masterstation'` (P3). If hourly values are available for another site in the vicinity of the station considered, the cumulative sub-daily mass curve can be transferred from the station that provides hourly values to the station of interest.

## Utilities
Among other features the `melodist.util` module includes some functions that might be useful for data analyses:

* `detect_gaps(dataframe, timestep, print_all=False, print_max=5, verbose=True)` can be used to find gaps in the data frame. A gap will be detected if any increment of time is not equal to the specified time step (in seconds).

* Some methods require time series with full days (24 h) only. `drop_incomplete_days(dataframe)` drops heading and tailing days if they are not complete (0-23h).

* For testing purposes an aggregation function is provided which aggregates hourly time series (data frames) to daily time series taking into account the characteristics of each meteorological variable (e.g., mean value for precipitation, daily total for precipitation, ...): `daily_from_hourly()`

## Data input/output

MELODIST includes a feature to read and save parameter settings for all disaggregation methods in order to transfer settings or to continue work at a later time. This feature is based on the JSON format which provides a flexible and easily readable ASCII file format for different applications.

To save MELODIST parameters included in station statistics object you can simply call the `to_json(filename)` method of this object. At any time it is possible to recall this settings by creating a new station statistics object based on the settings stored in that file:

```python
new_stationstatistics_object = melodist.StationStatistics.from_json(filename)
```

Since MELODIST is based on pandas, numerous ways to import and export pandas data frames exist. The `to_csv()` and `read_csv()` functions of pandas are ideal to load and save time series without any restriction with respect to MELODIST applications. 

MELODIST has also some additional specific data input/output capabilities in `melodist.data_io`, including functions to read data provided by the national weather services of the Netherlands (`read_single_knmi_file()`, `read_knmi_dataset()`) and Germany (`read_dwd()`). Moreover, the [SMET format](https://models.slf.ch/docserver/meteoio/SMET_specifications.pdf) is supported for reading and writing (`read_smet()`, `write_smet()`). This format is used in the [MeteoIO](https://models.slf.ch/p/meteoio) library.

## References
Förster, K., Hanzer, F., Winter, B., Marke, T., and Strasser, U.: An open-source MEteoroLOgical observation time series DISaggregation Tool (MELODIST v0.1.1), *Geosci. Model Dev.*, 9, 2315-2333, [doi:10.5194/gmd-9-2315-2016](https://doi.org/10.5194/gmd-9-2315-2016), 2016.

Hanzer, F., Förster, K., Nemec, J., and Strasser, U.: Projected cryospheric and hydrological impacts of 21st century climate change in the Ötztal Alps (Austria) simulated using a physically based approach, *Hydrol. Earth Syst. Sci.*, 22, 1593-1614, [doi:10.5194/hess-22-1593-2018](https://doi.org/10.5194/hess-22-1593-2018), 2018. 

## Version history

* 0.1.5 (03 Jul 2023):
  * Fix compatibility with recent numpy/pandas versions
* 0.1.4 (05 May 2020):
  * Upgraded to warrant pandas compatibility
* 0.1.3 (21 Nov 2018):
  * Sub-hourly precipitation disaggregation capabilities added along with a [new Jupyter notebook](examples/precip5min_example.ipynb)
  * bugfix in `precipitation.py`
  * speed-up of cascade statistics
* 0.1.2 (21 Dec 2017):
  * new disaggregation methods: `'sine_mean'`, `'mean_course_min_max'` and `'mean_course_mean'` for temperature, `'month_hour_precip_mean'` for humidity and `'mean_course'` for radiation
  * add option to preserve daily mean values in humidity disaggregation
  * add option to calculate and apply Angstroem/Bristow-Campbell parameters monthly or seasonally instead of for all values
  * bugfix in calculation of potential radiation
* 0.1.1 (03 Jun 2016):
  * data type corrections in order to avoid conversion warnings
  * new functions for the estimation of Angstrom and Bristow-Campbell parameters
  * handling of gaps in daily data fixed
  * speed-up of radiation computations
  * updates to example.py
  * pandas 0.18 API compatibility
* 0.1.0 (01 Mar 2016): First version of MELODIST
