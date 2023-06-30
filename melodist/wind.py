from __future__ import print_function, division, absolute_import
import melodist
import numpy as np
import pandas as pd
import scipy.optimize


def _cosine_function(x, a, b, t_shift):
    """genrates a diurnal course of windspeed accroding to the cosine function

    Args:
        x: series of euqally distributed windspeed values
        a: parameter a for the cosine function
        b: parameter b for the cosine function
        t_shift: parameter t_shift for the cosine function
        
    Returns:
        series including diurnal course of windspeed.
    """

    mean_wind, t = x
    return a * mean_wind * np.cos(np.pi * (t - t_shift) / 12) + b * mean_wind


def disaggregate_wind(wind_daily, method='equal', a=None, b=None, t_shift=None):
    """general function for windspeed disaggregation

    Args:
        wind_daily: daily values
        method: keyword specifying the disaggregation method to be used
        a: parameter a for the cosine function
        b: parameter b for the cosine function
        t_shift: parameter t_shift for the cosine function
        
    Returns:
        Disaggregated hourly values of windspeed.
    """
    assert method in ('equal', 'cosine', 'random'), 'Invalid method'

    wind_eq = melodist.distribute_equally(wind_daily)

    if method == 'equal':
        wind_disagg = wind_eq
    elif method == 'cosine':
        assert None not in (a, b, t_shift)
        wind_disagg = _cosine_function(np.array([wind_eq.values, wind_eq.index.hour]), a, b, t_shift)
    elif method == 'random':
        wind_disagg = wind_eq * (-np.log(np.random.rand(len(wind_eq))))**0.3

    return wind_disagg


def fit_cosine_function(wind):
    """fits a cosine function to observed hourly windspeed data

    Args:
        wind: observed hourly windspeed data
        
    Returns:
        parameters needed to generate diurnal features of windspeed using a cosine function
    """
    wind_daily = wind.groupby(wind.index.date).mean()
    wind_daily_hourly = pd.Series(index=wind.index, data=wind_daily.loc[wind.index.date].values)  # daily values evenly distributed over the hours

    df = pd.DataFrame(data=dict(daily=wind_daily_hourly, hourly=wind)).dropna(how='any')
    x = np.array([df.daily, df.index.hour])
    popt, pcov = scipy.optimize.curve_fit(_cosine_function, x, df.hourly)

    return popt
