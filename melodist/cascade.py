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
import numpy as np

"""
cascade.py provides basic definitions for the precipitation disaggregation
procedure using the cascade model by Olsson (1998)
"""

class BoxTypes:
    """BoxTypes represents the type of each box, which can be seen as one
    single value, with respect to its antecessor and follower along the time
    axis, respectively."""
    dry, starting, enclosed, ending, isolated = range(5)

class CascadeStatistics:
    """CascadeStatistics is a simple data structure including all relevant
    information based on statistical evaluation of hourly data in order to
    run the cascade based disaggregation applying these statistics for
    arbitrary daily precipitation time series."""
    def __init__(self):
        """constructor, initializes probability fields"""
        self.p01 = np.zeros((2,4))
        self.p10 = np.zeros((2,4))
        self.pxx = np.zeros((2,4))
        self.wxx = np.zeros((7,2,4))
        self.threshold = np.array([ 1.67093133,  2.46694444,  3.66730902,  5.39878419,  8.04924471])
        self.percentile = 50

    def fill_with_sample_data(self):
        """This function fills the corresponding object with sample data."""
        # replace these sample data with another dataset later
        # this function is deprecated as soon as a common file format for this
        # type of data will be available
        self.p01 = np.array([[0.576724636119866,  0.238722774405744,  0.166532122130638,  0.393474644666218],
                             [0.303345245644811,  0.0490956843857575, 0.0392403031072856, 0.228441890034704]])

        self.p10 = np.array([[0.158217002255554,  0.256581140990052,  0.557852226779526,  0.422638238585814],
                             [0.0439831163244427, 0.0474928027621488, 0.303675296728195,  0.217512052135178]])

        self.pxx = np.array([[0.265058361624580,  0.504696084604205,  0.275615651089836,  0.183887116747968],
                             [0.652671638030746,  0.903411512852094,  0.657084400164519,  0.554046057830118]])

        self.wxx = np.array([[[0.188389148850583, 0.0806836453984190, 0.0698113025807722, 0.0621499191745602],
                              [0.240993281622128, 0.0831019646519721, 0.0415130545715575, 0.155284541403192]],
                             [[0.190128959522795, 0.129220679033862,  0.0932213021787505, 0.193080698516532],
                              [0.196379692358065, 0.108549414860949,  0.0592714297292217, 0.0421945385836429]],
                             [[0.163043672107111, 0.152063537378127,  0.102823783410167,  0.0906028835221283],
                              [0.186579466868095, 0.189705690316132,  0.0990207345993082, 0.107831389238912]],
                             [[0.197765724699431, 0.220046257566978,  0.177876233348082,  0.261288786454262],
                              [0.123823472714948, 0.220514673922285,  0.102486496386323,  0.101975538893918]],
                             [[0.114435243444815, 0.170857634762767,  0.177327072603662,  0.135362730582518],
                              [0.0939211776723413,0.174291820501902,  0.125275822078525,  0.150842841725936]],
                             [[0.0988683809545079, 0.152323481100248, 0.185606883566286, 0.167242856061538],
                              [0.0760275616817939, 0.127275603247149,  0.202466168603738,  0.186580243138018]],
                             [[0.0473688704207573, 0.0948047647595988, 0.193333422312280,  0.0902721256884624],
                              [0.0822753470826286, 0.0965608324996108, 0.369966294031327,  0.255290907016382]]])

    def __add__(self, other_casc_obj):
        self.p01 = self.p01[:,:] + other_casc_obj.p01[:,:]
        self.p10 = self.p10[:,:] + other_casc_obj.p10[:,:]
        self.pxx = self.pxx[:,:] + other_casc_obj.pxx[:,:]
        for k in range(0,7):
            self.wxx[k,:,:] = self.wxx[k,:,:] + other_casc_obj.wxx[k,:,:]

    def __mul__(self, factor):
        self.p01 = self.p01[:,:] * factor
        self.p10 = self.p10[:,:] * factor
        self.pxx = self.pxx[:,:] * factor
        self.wxx = self.wxx[:,:,:] * factor

    @classmethod
    def from_dict(cls, d):
        casc = cls()

        casc.p01 = np.array(d['p01'])
        casc.p10 = np.array(d['p10'])
        casc.pxx = np.array(d['pxx'])
        casc.wxx = np.array(d['wxx'])
        casc.threshold = np.array(d['threshold'])
        casc.percentile = d['percentile']

        assert casc.p01.shape == (2, 4)
        assert casc.p10.shape == (2, 4)
        assert casc.pxx.shape == (2, 4)
        assert casc.wxx.shape == (7, 2, 4)

        return casc
