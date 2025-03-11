#!/usr/bin/env python

#This file is part of chelsa_cmip6.
#
#chelsa_cmip6 is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.

#chelsa_cmip6 is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.

#You should have received a copy of the GNU General Public License
#along with chelsa_cmip6.  If not, see <https://www.gnu.org/licenses/>.


import xarray as xr
import numpy as np
from scipy.interpolate import interp1d


def growing_degree_days(tas, threshold=None):
    """ 
    Calculate growing degree days
    
    :param tas: Daily mean 2m air temperature [Celsius]
    :param threshold: The threshold temperature above which a day is counted as growing day [Celsius]

    :return: growing degree days [Celsius]
    :rtype: float
    """

    if threshold == None:
        threshold = 5

    if len(tas) == 366 or len(tas) == 365:
        gdd = np.sum([i for i in tas if i >= threshold])

    if len(tas) == 12:
        # JAN FEB MAR APR MAY JUN JUL AUG SEP OCT NOV DEC
        # 0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334
        midmonth = [-15, 15, 45, 74, 105, 135, 166, 196, 227, 258, 288, 319, 349, 380]
        monthv = [12, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 1]

        tas14 = []
        for n in range(0, 14):
            tas14.append(tas[monthv[n] - 1])

        f2 = interp1d(midmonth, tas14, kind='cubic')
        xnew = np.linspace(0, 365, num=366, endpoint=True)
        tas365 = f2(xnew)
        gdd = np.sum([i for i in tas365 if i >= threshold])

    return gdd


class quarter_class:
    """ 
    quarters class for monthly climatologies 
    
    Creates quarters from monthly climatologies.

    :param target_variable: the variable that is to be aggregated for a quarter_class
    :param quarter_variable: the variable that forms the basis for the quarters (e.g. precipitation for 'wettest' quarter)
    :param agg_target: the way how the target_variable is aggregated (e.g. mean)
    :param agg_quarter: the way how the quarter_variable is aggregated (e.g. max)
    :param find_fun: determines what needs to be found (e.g. minimum temperature of coldes quarter = 'min')

    """
    def __init__(self, target_variable, quarter_variable, agg_target, agg_quarter, find_fun):
        self.target_variable = target_variable
        self.quarter_variable = quarter_variable
        self.agg_target = agg_target #how should array1 be aggregated
        self.agg_quarter = agg_quarter #how should array2 be aggregated
        self.find_fun = find_fun #is the min or the max whats need to be found

    #[dec, jan, feb, mar, apr, may, jun, jul, sep, oct, nov, dec, jan]
    def _create_quarter_(self, xv, agg):
        """
        creates quarters from a vector of 12 monthly means

        :param xv: monthly values of a variable
        :param agg: aggregation method

        :return: array of size 4 containing the quarters
        :rtype: array
        """
        #xv = [11, 20, 30, 104, 95, 96, 75, 85, 90, 190, 181, 172]
        monthv = [12, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 1]
        b1 = []
        for n in range(0, 14):
            b1.append(xv[monthv[n]-1])
        a1 = []
        for m in range(1, 13):
            global x0
            if agg == 'sum':
                x0 = np.sum([b1[m - 1], b1[m], b1[m + 1]])
            if agg == 'mean':
                x0 = np.mean([b1[m - 1], b1[m], b1[m + 1]])
            if agg == 'max':
                x0 = np.max([b1[m - 1], b1[m], b1[m + 1]])
            if agg == 'min':
                x0 = np.min([b1[m - 1], b1[m], b1[m + 1]])
            a1.append(x0)
        return a1

    def _comp_quarters_(self, target_variable, quarter_variable): # ,array1, array2, agg1, agg2, fun2
        """
        Compare quarters for a target variable

        :param target_variable: the variable that is to be aggregated for a quarter_class
        :param quarter_variable: the variable that forms the basis for the quarters (e.g. precipitation for 'wettest' quarter)

        :return: positional location of the comparision (e.g. the position of the coldest quarter in the array)
        :rtype: integer
        """
        target_quarter = self._create_quarter_(target_variable, self.agg_target)
        quarter_quarter = self._create_quarter_(quarter_variable, self.agg_quarter)
        if self.find_fun == 'max':
            a1 = np.max(quarter_quarter)
        if self.find_fun == 'min':
            a1 = np.min(quarter_quarter)
        a = target_quarter[quarter_quarter.index(a1)]
        return a

    def comp_quarters_array(self):
        """
        Compare quarters over an xarray

        :return: array of target variable for a quarter
        :rtype: xarray
        """
        res_arr = xr.apply_ufunc(self._comp_quarters_,  # function to apply
                                 self.target_variable,
                                 self.quarter_variable,# pass arguments.
                                 input_core_dims=[['month'], ['month']],
                                 vectorize=True,
                                 dask='parallelized',  # let dask handle the parallelization
                                 dask_gufunc_kwargs=['allow_rechunk'],
                                 output_dtypes=[np.float32])  # data type of the output(s)
        return res_arr


class BioClim:
    """ 
    Climatology class to calculate bioclimatic variables from monthly climatologies of temperatures
    and precipitation. 

    :param tas: Monthly mean daily mean 2m air temperature [Celsius]
    :param tasmax: Monthly mean daily maximum 2m air temperature [Celsius]
    :param tasmin: Monthly mean daily minimum 2m air temperature [Celsius]
    :param pr: Monthly precipitation rate [kg*m**-2*month**-1]
    """
    def __init__(self, pr, tas, tasmax, tasmin):
        """ Create a set of baseline clims """
        self.tas = tas.load() #chunk({month: -1}) #rename({'tas': 'var'})
        self.tasmax = tasmax.load() #.chunk({month: -1}) #rename({'tasmax': 'var'})
        self.tasmin = tasmin.load() #.chunk({month: -1}) #rename({'tasmin': 'var'})
        self.pr = pr.load() #chunk({month: -1}) #rename({'pr': 'var'})

    def _mean_(self, x ):
        """
        calculate mean

        :param x: a variable 
        :return: mean
        :rtype: float
        """
        s = np.sum(x)
        n = len(x)
        mean = s/n
        return mean

    def _diurnalrange_(self, tasmax, tasmin):
        """
        Calculate the mean monthly diurnal range of temperatures
        
        :param tasmax: Monthly mean daily maximum 2m air temperature [Celsius]
        :param tasmin: Monthly mean daily minimum 2m air temperature [Celsius]
        :return: diurnal range of temperatures [Kelvin]
        :rtype: float
        """
        return np.sum(tasmax - tasmin) / 12

    def _sd_(self, x):
        """
        calculate standard deviation

        :param x: a variable 
        :return: standard deviation
        :rtype: float
        """
        return np.std(x)

    def _max_(self, x):
        """
        calculate maximum

        :param x: a variable 
        :return: maximum
        :rtype: float
        """
        return np.max(x)

    def _min_(self, x):
        """
        calculate minimum

        :param x: a variable 
        :return: minimum
        :rtype: float
        """
        return np.min(x)

    def _sum_(self, x):
        """
        calculate sum

        :param x: a variable 
        :return: sum
        :rtype: float
        """
        return np.sum(x)

    def _cv_(self, x):
        """
        calculate coefficient or variation

        :param x: a variable 
        :return: coefficient of variation
        :rtype: float
        """
        sigma = self._sd_(x)
        mu = self._mean_(x)
        cv = sigma*100/mu
        return cv

    def _bio3_(self, tasmax, tasmin):
        """
        calculate bio3 "Isothermality"
        The ratio of diurnal variation (bio2) to annual range of air temperature (bio7)

        :param tasmax: Monthly mean daily maximum 2m air temperature [Celsius]
        :param tasmin: Monthly mean daily minimum 2m air temperature [Celsius]
        :return: Isothermality [unitless]
        :rtype: float
        """
        bio5 = self._max_(tasmax)
        bio6 = self._min_(tasmin)
        bio7 = bio5-bio6
        bio2 = self._diurnalrange_(tasmax, tasmin)
        bio3 = bio2 / bio7
        return bio3

    def _bio7_(self, tasmax, tasmin):
        """
        calculate bio7 "annual range of air temperature"
        The difference between the Maximum Temperature of Warmest month and the Minimum Temperature of Coldest month

        :param tasmax: Monthly mean daily maximum 2m air temperature [Celsius]
        :param tasmin: Monthly mean daily minimum 2m air temperature [Celsius]
        :return: annual range of air temperature [Kelvin]
        :rtype: float
        """
        bio5 = self._max_(tasmax)
        bio6 = self._min_(tasmin)
        bio7 = bio5-bio6
        return bio7

    def bio1(self):
        """
        Create mean annual temperature

        :return: mean annual 2m air temperature [Celsius]
        :rtype: xarray
        """
        res_arr = xr.apply_ufunc(self._mean_,  # function to apply
                                 self.tas['tas'],  # pass arguments.
                                 input_core_dims=[['month']],
                                 vectorize=True,
                                 dask='parallelized',
                                 dask_gufunc_kwargs=['allow_rechunk'],
                                 output_dtypes=[np.float32]) #dask='parallelized',  # let dask handle the parallelization
                                  # data type of the output(s)
        res_arr = res_arr.to_dataset(name='bio1')
        return res_arr

    def bio2(self):
        """
        Create mean diurnal temperature range

        :return: mean diurnal temperature range [Kelvin]
        :rtype: xarray
        """
        res_arr = xr.apply_ufunc(self._diurnalrange_,
                                 self.tasmax['tasmax'], self.tasmin['tasmin'],
                                 input_core_dims=[['month'], ['month']],
                                 vectorize=True,
                                 dask='parallelized',
                                 dask_gufunc_kwargs=['allow_rechunk'],
                                 output_dtypes=[np.float32])
        res_arr = res_arr.to_dataset(name='bio2')
        return res_arr

    def bio3(self):
        """ 
        Isothermality

        :return: Isothermality [unitless]
        :rtype: xarray
        """
        res_arr = xr.apply_ufunc(self._bio3_,
                                 self.tasmax['tasmax'],
                                 self.tasmin['tasmin'],
                                 input_core_dims=[['month'],['month']],
                                 dask='parallelized',
                                 dask_gufunc_kwargs=['allow_rechunk'],
                                 vectorize=True,
                                 output_dtypes=[np.float32])
        res_arr = res_arr.to_dataset(name='bio3')
        return res_arr

    def bio4(self):
        """
        Temperature Seasonality (Standard Deviation)

        :return: Temperature Seasonality
        :rtype: xarray
        """
        res_arr = xr.apply_ufunc(self._sd_,
                                 self.tas['tas'],
                                 input_core_dims=[['month']],
                                 vectorize=True,
                                 dask='parallelized',
                                 dask_gufunc_kwargs=['allow_rechunk'],
                                 output_dtypes=[np.float32])
        res_arr = res_arr.to_dataset(name='bio4')
        return res_arr

    def bio5(self):
        """
        Max Temperature of Warmest Month

        :return: Maximum temperature of the warmest month [Celsius]
        :rtype: xarray
        """
        res_arr = xr.apply_ufunc(self._max_,
                                 self.tasmax['tasmax'],
                                 input_core_dims=[['month']],
                                 vectorize=True,
                                 dask='parallelized',
                                 dask_gufunc_kwargs=['allow_rechunk'],
                                 output_dtypes=[np.float32])
        res_arr = res_arr.to_dataset(name='bio5')
        return res_arr

    def bio6(self):
        """
        Minimum Temperature of Coldest Month

        :return: Minimum temperature of the coldest month [Celsius]
        :rtype: xarray
        """
        res_arr = xr.apply_ufunc(self._min_,
                                 self.tasmin['tasmin'],
                                 input_core_dims=[['month']],
                                 vectorize=True,
                                 dask='parallelized',
                                 dask_gufunc_kwargs=['allow_rechunk'],
                                 output_dtypes=[np.float32])
        res_arr = res_arr.to_dataset(name='bio6')
        return res_arr

    def bio7(self):
        """
        Annual Temperature Range
        
        :return: Annual range in temperature [Celsius]
        :rtype: xarray
        """
        res_arr = xr.apply_ufunc(self._bio7_,
                                 self.tasmax['tasmax'], self.tasmin['tasmin'],
                                 input_core_dims=[['month'], ['month']],
                                 vectorize=True,
                                 dask='parallelized',
                                 dask_gufunc_kwargs=['allow_rechunk'],
                                 output_dtypes=[np.float32])
        res_arr = res_arr.to_dataset(name='bio7')
        return res_arr

    def bio8(self):
        """
        Mean Temperature of Wettest Quarter
        
        :return: Mean temperature of the wettest quarter [Celsius]
        :rtype: xarray
        """
        res_arr = quarter_class(target_variable=self.tas['tas'],
                                quarter_variable=self.pr['pr'],
                                agg_target="mean",
                                agg_quarter="sum",
                                find_fun="max").comp_quarters_array()
        res_arr = res_arr.to_dataset(name='bio8')
        return res_arr

    def bio9(self):
        """
        Mean Temperature of the driest Quarter
        
        :return: Mean temperature of the driest quarter [Celsius]
        :rtype: xarray
        """
        res_arr = quarter_class(target_variable=self.tas['tas'],
                                quarter_variable=self.pr['pr'],
                                agg_target="mean",
                                agg_quarter="sum",
                                find_fun="min").comp_quarters_array()
        res_arr = res_arr.to_dataset(name='bio9')
        return res_arr

    def bio10(self):
        """
        Mean Temperature of the warmest Quarter
        
        :return: Mean temperature of the warmest quarter [Celsius]
        :rtype: xarray
        """
        res_arr = quarter_class(target_variable=self.tas['tas'],
                                quarter_variable=self.tas['tas'],
                                agg_target="mean",
                                agg_quarter="mean",
                                find_fun="max").comp_quarters_array()
        res_arr = res_arr.to_dataset(name='bio10')
        return res_arr

    def bio11(self):
        """
        Mean Temperature of coldest Quarter
        
        :return: Mean temperature of the coldest quarter [Celsius]
        :rtype: xarray
        """
        res_arr = quarter_class(target_variable=self.tas['tas'],
                                quarter_variable=self.tas['tas'],
                                agg_target="mean",
                                agg_quarter="mean",
                                find_fun="min").comp_quarters_array()
        res_arr = res_arr.to_dataset(name='bio11')
        return res_arr

    def bio12(self):
        """
        Annual Precipitation Sum
        
        :return: Annual precipition sum [kg*m**-2*year*-1]
        :rtype: xarray
        """
        res_arr = xr.apply_ufunc(self._sum_,
                                 self.pr['pr'],
                                 input_core_dims=[['month']],
                                 vectorize=True,
                                 dask='parallelized',
                                 dask_gufunc_kwargs=['allow_rechunk'],
                                 output_dtypes=[np.float32])
        res_arr = res_arr.to_dataset(name='bio12')
        return res_arr

    def bio13(self):
        """Precipitation of wettest month"""
        res_arr = xr.apply_ufunc(self._max_,
                                 self.pr['pr'],
                                 input_core_dims=[['month']],
                                 vectorize=True,
                                 dask='parallelized',
                                 dask_gufunc_kwargs=['allow_rechunk'],
                                 output_dtypes=[np.float32])
        res_arr = res_arr.to_dataset(name='bio13')
        return res_arr

    def bio14(self):
        """
        Precipitation of driest month
                
        :return: precipition rate [kg*m**-2*month*-1]
        :rtype: xarray
        """
        res_arr = xr.apply_ufunc(self._min_,
                                 self.pr['pr'],
                                 input_core_dims=[['month']],
                                 vectorize=True,
                                 dask='parallelized',
                                 dask_gufunc_kwargs=['allow_rechunk'],
                                 output_dtypes=[np.float32])
        res_arr = res_arr.to_dataset(name='bio14')
        return res_arr

    def bio15(self):
        """
        Precipitation Seasonality  

        :return: coefficient of variation in precipition rate [kg*m**-2*month*-1]
        :rtype: xarray
        """
        res_arr = xr.apply_ufunc(self._cv_,
                                 self.pr['pr'],
                                 input_core_dims=[['month']],
                                 vectorize=True,
                                 dask='parallelized',
                                 dask_gufunc_kwargs=['allow_rechunk'],
                                 output_dtypes=[np.float32])
        res_arr = res_arr.to_dataset(name='bio15')
        return res_arr

    def bio16(self):
        """
        Precipitation of Wettest Quarter    

        :return: precipition rate [kg*m**-2*month*-1]
        :rtype: xarray
        """
        res_arr = quarter_class(target_variable=self.pr['pr'],
                                quarter_variable=self.pr['pr'],
                                agg_target="sum",
                                agg_quarter="sum",
                                find_fun="max").comp_quarters_array()
        res_arr = res_arr.to_dataset(name='bio16')
        return res_arr

    def bio17(self):
        """
        Precipitation of Driest Quarter                
        
        :return: precipition rate [kg*m**-2*month*-1]
        :rtype: xarray
        """
        res_arr = quarter_class(target_variable=self.pr['pr'],
                                quarter_variable=self.pr['pr'],
                                agg_target="sum",
                                agg_quarter="sum",
                                find_fun="min").comp_quarters_array()
        res_arr = res_arr.to_dataset(name='bio17')
        return res_arr

    def bio18(self):
        """
        Precipitation of Warmest Quarter 

        :return: precipition rate [kg*m**-2*month*-1]
        :rtype: xarray
        """
        res_arr = quarter_class(target_variable=self.pr['pr'],
                                quarter_variable=self.tasmax['tasmax'],
                                agg_target="sum",
                                agg_quarter="mean",
                                find_fun="max").comp_quarters_array()
        res_arr = res_arr.to_dataset(name='bio18')
        return res_arr

    def bio19(self):
        """
        Precipitation of Coldest Quarter                
        
        :return: precipition rate [kg*m**-2*month*-1]
        :rtype: xarray
        """
        res_arr = quarter_class(target_variable=self.pr['pr'],
                                quarter_variable=self.tas['tas'],
                                agg_target="sum",
                                agg_quarter="mean",
                                find_fun="min").comp_quarters_array()
        res_arr = res_arr.to_dataset(name='bio19')
        return res_arr

    def gdd(self):
        """
        Growing degree days
                    
        :return: groeing degree days [Kelvin]
        :rtype: xarray
        """
        res_arr = xr.apply_ufunc(growing_degree_days,
                                 self.tas['tas'],
                                 input_core_dims=[['month']],
                                 vectorize=True,
                                 dask='parallelized',
                                 dask_gufunc_kwargs=['allow_rechunk'],
                                 output_dtypes=[np.float32])
        res_arr = res_arr.to_dataset(name='gdd')
        return res_arr

