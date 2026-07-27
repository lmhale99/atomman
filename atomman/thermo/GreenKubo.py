from typing import Optional, Union, Tuple

import numpy as np
import numpy.typing as npt

import scipy

# Local imports
import atomman.unitconvert as uc

############################ Supporting functions ############################

def rolling_mean(x: npt.ArrayLike,
                 delta: int):
    """
    Evaluates the rolling mean of a series.

    Parameters
    ----------
    x : array-like
        The 1D data series to evaluate.
    delta : int
        The rolling window size for evaluating the mean values.

    Returns
    -------
    mean : numpy.ndarray
        An array containing the mean values.  Will have (delta-1) less values
        than x.
    """
    # Create a sliding view of the data
    slide_view = np.lib.stride_tricks.sliding_window_view(x, delta)
    
    # Compute the mean of each window
    mean = np.mean(slide_view, axis=-1)

    return mean

def rolling_std(x: npt.ArrayLike,
                delta: int,
                batchsize: int = 100):
    """
    Evaluates the rolling standard deviation of a series.

    Parameters
    ----------
    x : array-like
        The 1D data series to evaluate.
    delta : int
        The rolling window size for evaluating the std values.
    batchsize : int, optional
        Indicates how many std evaluations to make at a time.  The batchsize
        only affects memory usage and computational time.  I think numpy.std()
        creates a temporary array in memory which can be large if the size of x
        is large.  Keeping the batchsize small avoids memory issues and speeds
        up execution.  Default value is 100.

    Returns
    -------
    std : numpy.ndarray
        An array containing the std values.  Will have (delta-1) less values
        than x.
    """
    # Create a sliding view of the data
    slide_view = np.lib.stride_tricks.sliding_window_view(x, delta)

    # Batch compute the std of each window
    std = np.empty(slide_view.shape[0])
    nbatches = std.shape[0] // batchsize + 1
    for i in range(nbatches):
        s = i * batchsize
        e = s + batchsize
        std[s:e] = np.std(slide_view[s:e], axis=-1)

    # Return the rolling std
    return std

def autocorrelationfunction(x: npt.ArrayLike):
    """
    Computes the auto-correlation function for a 1D sequence.

    Parameters
    ----------
    x : array-like
        The 1D sequence to evaluate the auto-correlation function for.
    """
    
    x = np.asarray(x)
    assert x.ndim == 1, 'input must be a 1D array'
    
    # Compute the full auto-correlation sums
    xx = scipy.signal.correlate(x, x, mode='full')[-len(x):]

    # Divide by count to get average
    xx /= np.arange(x.shape[0], 0, -1)
    
    return xx

def cumtrapezoid(y: npt.ArrayLike,
                 x: Optional[npt.ArrayLike] = None,
                 dx: Union[float, npt.ArrayLike] = 1.0):
    """
    Computes the cumulative trapezoid integral of a series.

    Parameters
    ----------
    y : array-like
        The 1D data series.
    x : array-like, optional
        The sample points corresponding to the y values. If x is None, the
        sample points are assumed to be evenly spaced dx apart. The default is
        None.
    dx : float, optional
        The corresponding x spacing between all y values.  Default value is 1,
        i.e. no scaling.
    """
    # Create dx array if x is given
    if x is not None:
        x = np.asarray(x, dtype=float)
        dx = x[1:] - x[:-1]
    
    # Create empty output array and set first value
    y = np.asarray(y, dtype=float)
    cumintegral = np.empty_like(y, dtype=float)
    cumintegral[0] = 0.0

    # Compute individual trapezoid areas and cumulatively sum it
    cumintegral[1:] = np.cumsum(dx * (y[1:] + y[:-1]) / 2)
    
    return cumintegral


######################### Green-Kubo manager classes #########################

class GreenKubo():

    def __init__(self,
                 time: npt.ArrayLike,
                 x: Optional[npt.ArrayLike] = None,
                 acf: Optional[npt.ArrayLike] = None,
                 acfsamples: Union[int, npt.ArrayLike] = 1):
        """
        Base class for Green-Kubo calculations.  This takes the core
        input data series, computes the auto-correlation, fluctuation,
        and the cumulative integral for the data.
        
        Parameters
        ----------
        time : array-like
            The values of time for the 1D array or acf.
        x : array-like or None, optional
            The 1D array that the auto-correlation function is computed for, i.e.
            <x(t)*x(0)>. Either value or acf must be given.
        acf : array-like or None, optional
            A pre-computed auto-correlation function to evaluate.  Either
            x or acf must be given.
        acfsamples : int or array-like, optional
            If acf is given this indicates how many (relative) samples were
            averaged over for each given acf value.  This is typically either a
            constant (default = 1) if the same number of samples were made for each
            acf, or is a linearly decreasing array if the acf was computed for an
            ongoing known time series value.  This is automatically determined if
            x is given.
        """
        self.__time = np.asarray(time, dtype=float)
        if self.time.ndim != 1:
            raise ValueError('time must be a 1D array')

        if x is not None:
            if acf is not None:
                raise ValueError('cannot give x and acf!')
            
            # Set and check x values
            self.__x = np.asarray(x, dtype=float)
            self.__acfsamples = np.arange(self.x.shape[0], 0, -1)

            if self.x.ndim != 1:
                raise ValueError('x must be a 1D array')
            if self.time.shape[0] != self.x.shape[0]:
                raise ValueError(f'time and x have different number of values: {self.time.shape[0]} != {self.x.shape[0]}')

            # Compute the auto-correlation function of x
            self.__acf = autocorrelationfunction(self.x)
            
        else:
            # Set and check acf values
            self.__x = None
            self.__acf = np.asarray(acf, dtype=float)
            if isinstance(acfsamples, int):
                self.__acfsamples = np.full(self.acf.shape, acfsamples)
            else:
                self.__acfsamples = np.asarray(acfsamples)
            if self.acf.ndim != 1:
                raise ValueError('acf must be a 1D array')
            if self.time.shape[0] != self.acf.shape[0]:
                raise ValueError(f'time and acf have different number of values: {self.time.shape[0]} != {self.acf.shape[0]}')

        # Compute the cumulative integral
        self.__integral = cumtrapezoid(self.acf, self.time)

    @classmethod
    def acf_from_lammps(cls,
                        filename: str,
                        name: str,
                        timestep: float,
                        unit: Optional[str] = None,
                        index: int = 2147483647):
        """
        Read in step and auto-correlation data from LAMMPS fix ave/correlate
        output.

        Parameters
        ----------
        filename : str or Path
            The data file generated by the fix ave/correlate file.
        name : str
            The column name associated with the auto-correlation function data.
        timestep : float
            The timestep used by the LAMMPS simulation, i.e.
            time is TimeDelta * timestep.
        unit : str or None, optional
            The units that the acf data is in. If None (default), no unit
            conversions will be done.
        index : int, optional
            The fix ave/correlate file may have multiple time series runs.
            This allows for the selection of which series to read in.  If
            index is greater than the number of time series runs in the file,
            the final series will be read in (which is what is usually wanted).
            The default value is 2147483647, the 32-bit maximum int, which
            should in a practical sense always be larger than the number of series.
        """
        # Read the file
        with open(filename) as f:
            lines = f.readlines()
        
        # Read the column header
        terms = lines[2].split()
        step_index = terms.index('TimeDelta') - 1
        data_index = terms.index(name) - 1
       
        i = 3
        while len(lines) > i:
            terms = lines[i].split()
            if len(terms) == 0:
                break
            count = int(terms[1])
            i += count + 1
            index -= 1
            if index < 0:
                break

        time = np.empty(count) 
        acf = np.empty(count)

        for j, line in enumerate(lines[i-count:i]):
            terms = line.split()
            time[j] = float(terms[step_index])
            acf[j] = float(terms[data_index])

        # Convert time and acf values
        time = time * timestep
        acf = uc.set_in_units(acf, unit)

        return cls(time=time, acf=acf)

    @property
    def x(self) -> Optional[np.ndarray]:
        """numpy.ndarray or None: The input array for the auto-correlation function"""
        return self.__x

    @property
    def time(self) -> np.ndarray:
        """numpy.ndarray: The time values"""
        return self.__time

    @property
    def acf(self) -> np.ndarray:
        """numpy.ndarray: The auto-correlation function"""
        return self.__acf
    
    @property
    def acfsamples(self) -> np.ndarray:
        """numpy.ndarray: The number of samples that each acf value was averaged over"""
        return self.__acfsamples
        
    @property
    def integral(self) -> np.ndarray:
        """numpy.ndarray: The cumulative integral of the auto-correlation function"""
        return self.__integral

    def fluctuation(self,
                    delta: int,
                    batchsize: int = 100) -> np.ndarray:
        """
        Evaluates the "fluctuation" function of the auto-correlation function
        evaluating |std(acf) / mean(acf)| for rolling windows.  As the data's
        signal is related to its mean and the noise is related to its std, this
        can be used to identify where the noise in the data is larger than the
        signal.

        Parameters
        ----------
        delta : int
            The rolling window size for evaluating the fluctuation.
        batchsize : int, optional
            Indicates how many std evaluations to make at a time.  The batchsize
            only affects memory usage and computational time.  I think numpy.std()
            creates a temporary array in memory which can be large if the size of x
            is large.  Keeping the batchsize small avoids memory issues and speeds
            up execution.  Default value is 100.

        Returns
        -------
        fluctuation : numpy.ndarray
            The fluctuation of the acf.  This array will have (delta-1) fewer
            values than the acf array.
        """
        mean = rolling_mean(self.acf, delta)
        std = rolling_std(self.acf, delta, batchsize=batchsize)

        # Return the rolling |std/mean| values
        return np.abs(std / mean)

    def std_noise(self,
                  startindex: Optional[int] = None,
                  starttime: Optional[float] = None):
        """
        Evaluates the standard deviation of the tail end of the acf data to
        obtain an estimate of the "pure noise".

        Parameters
        ----------
        startindex : int or None, optional
            All acf data starting with this index will be used to compute the
            standard deviation.  This value should be sufficiently large such
            that the acf at that point is pure noise.  Cannot be given with
            starttime.  If neither startindex or starttime are given,
            startindex will be taken as the midpoint of the acf array.
        starttime : float or None, optional
            All acf data starting with this time value will be used to compute
            the standard deviation.  This value should be sufficiently large such
            that the acf at that point is pure noise.  Cannot be given with
            startindex.  If neither startindex or starttime are given,
            startindex will be taken as the midpoint of the acf array.

        Returns
        -------
        std : numpy.ndarray
            The time-dependent standard deviation for the acf data.  If acfsamples
            is constant, then std is constant.  Otherwise, std is scaled by the
            number of acf samples.
        """
        # Handle startindex/starttime values
        if startindex is not None:
            if starttime is not None:
                raise ValueError('startindex and starttime cannot both be given')
        elif starttime is not None:
            startindex = np.where(self.time >= starttime)[0][0]
        else:
            startindex = self.acf.shape[0] // 2

        # Normalize acf values using the number of samples
        sqrtN = self.acfsamples**0.5
        const = np.std(self.acf[startindex:] * sqrtN[startindex:]**0.5)

        # Scale the std constant using the number of samples
        return const / sqrtN

    def std_noise_fluctuation(self,
                              delta: int,
                              startindex: Optional[int] = None,
                              starttime: Optional[float] = None) -> np.ndarray:
        """
        Evaluates the "fluctuation" function of the auto-correlation function
        evaluating |std(acf) / mean(acf)| using a rolling window for the mean
        and the "pure noise" standard deviation.  As the data's signal is
        related to its mean and the noise is related to its std, this
        can be used to identify where the noise in the data is larger than the
        signal.

        Parameters
        ----------
        delta : int
            The rolling window size for evaluating the fluctuation.
        startindex : int or None, optional
            All acf data starting with this index will be used to compute the
            standard deviation.  This value should be sufficiently large such
            that the acf at that point is pure noise.  Cannot be given with
            starttime.  If neither startindex or starttime are given,
            startindex will be taken as the midpoint of the acf array.
        starttime : float or None, optional
            All acf data starting with this time value will be used to compute
            the standard deviation.  This value should be sufficiently large such
            that the acf at that point is pure noise.  Cannot be given with
            startindex.  If neither startindex or starttime are given,
            startindex will be taken as the midpoint of the acf array.

        Returns
        -------
        fluctuation : numpy.ndarray
            The fluctuation of the acf.  This array will have (delta-1) fewer
            values than the acf array.
        """
        mean = rolling_mean(self.acf, delta)

        # Compute the std of the pure noise region
        std = self.std_noise(startindex=startindex, starttime=starttime)

        # Trim std to the same length as mean
        std = std[:mean.shape[0]]

        # Return the rolling |std/mean| values
        return np.abs(std / mean)
    
    def std_noise_fraction(self,
                           delta: int,
                           startindex: Optional[int] = None,
                           starttime: Optional[float] = None) -> np.ndarray:
        """
        Computes the fraction of acf values within a given rolling window
        that are less than the "pure noise" standard deviation.  Large
        fractions indicate that most acf values are indistinguishable
        from the noise.

        Parameters
        ----------
        delta : int
            The size of the rolling window to use.  
        startindex : int or None, optional
            All acf data starting with this index will be used to compute the
            standard deviation.  This value should be sufficiently large such
            that the acf at that point is pure noise.  Cannot be given with
            starttime.  If neither startindex or starttime are given,
            startindex will be taken as the midpoint of the acf array.
        starttime : float or None, optional
            All acf data starting with this time value will be used to compute
            the standard deviation.  This value should be sufficiently large such
            that the acf at that point is pure noise.  Cannot be given with
            startindex.  If neither startindex or starttime are given,
            startindex will be taken as the midpoint of the acf array.
        
        Returns
        -------
        fraction : numpy.ndarray
            The fractional count of acf points less than the pure noise std
            for each window.  This array will have (delta-1) fewer
            values than the acf array.
        """
        # Compute the std of the pure noise region
        std = self.std_noise(startindex=startindex, starttime=starttime)

        # Compare |acf| to std for all points
        withinstd = np.abs(self.acf) < std

        # Rolling mean of booleans gives fractional count of Trues in each window
        fraction = rolling_mean(withinstd, delta)

        return fraction

    def tcut_fluctuation(self,
                         delta: int,
                         threshold: Optional[float] = None,
                         batchsize: int = 100,
                         timeshift: str = 'first') -> Tuple[int, float]:
        """
        Identifies a cutoff value for the calculation based on the fluctuation
        method.

        Parameters
        ----------
        delta : int
            The rolling window size for evaluating the fluctuation.
        threshold : float or None, optional
            The threshold value to use.  The cutoff is identified based on the
            first acf value where the fluctuation is greater than the threshold
            value.  Default value of None will set the threshold to the max of 1
            or twice the fluctuation value at t=0. 
        batchsize : int, optional
            Indicates how many std evaluations to make at a time.  The batchsize
            only affects memory usage and computational time.  I think numpy.std()
            creates a temporary array in memory which can be large if the size of x
            is large.  Keeping the batchsize small avoids memory issues and speeds
            up execution.  Default value is 100.
        timeshift : str, optional
            Indicates where the time value returned is taken from within the
            rolling window: "first" (default), "middle" or "last".

        Returns
        -------
        icut : int
            The array index associated with the cutoff value.
        tcut : float
            The time value associated with the cutoff value.
        """
        # Compute the fluctuation data
        f = self.fluctuation(delta, batchsize=batchsize)

        # Set the default threshold
        if threshold is None:
            threshold = max(1, 2*f[0])

        # Identify the first index where f > threshold
        icut: int = np.where(f >= threshold)[0][0]

        # Shift according to the timeshift option
        if timeshift == 'last':
            icut += delta - 1
        elif timeshift == 'middle':
            icut += (delta - 1) // 2
        elif timeshift != 'first':
            raise ValueError('invalid timeshift option')
        
        return icut, self.time[icut]

    def tcut_std_noise_fluctuation(self,
                                   delta: int,
                                   threshold: float = 1.0,
                                   startindex: Optional[int] = None,
                                   starttime: Optional[float] = None,
                                   timeshift: str = 'first') -> Tuple[int, float]:
        """
        Identifies a cutoff value for the calculation based on the fluctuation
        method and the pure noise std.

        Parameters
        ----------
        delta : int
            The size of the rolling window to use.  
        threshold : float or None, optional
            The threshold value to use.  The cutoff is identified based on the
            first acf value where the fluctuation is greater than the threshold
            value.  Default value is 1.0 (i.e. std == mean). 
        startindex : int or None, optional
            All acf data starting with this index will be used to compute the
            standard deviation.  This value should be sufficiently large such
            that the acf at that point is pure noise.  Cannot be given with
            starttime.  If neither startindex or starttime are given,
            startindex will be taken as the midpoint of the acf array.
        starttime : float or None, optional
            All acf data starting with this time value will be used to compute
            the standard deviation.  This value should be sufficiently large such
            that the acf at that point is pure noise.  Cannot be given with
            startindex.  If neither startindex or starttime are given,
            startindex will be taken as the midpoint of the acf array.
        timeshift : str, optional
            Indicates where the time value returned is taken from within the
            rolling window: "first" (default), "middle" or "last".
        """
        # Compute the fluctuation data
        f = self.std_noise_fluctuation(delta, startindex=startindex, starttime=starttime)

        # Set the default threshold
        if threshold is None:
            threshold = max(1, 2*f[0])

        # Identify the first index where f > threshold
        icut: int = np.where(f >= threshold)[0][0]

        # Shift according to the timeshift option
        if timeshift == 'last':
            icut += delta - 1
        elif timeshift == 'middle':
            icut += (delta - 1) // 2
        elif timeshift != 'first':
            raise ValueError('invalid timeshift option')
        
        return icut, self.time[icut]

    def tcut_std_noise_fraction(self,
                                delta: int,
                                threshold: float = 0.8,
                                startindex: Optional[int] = None,
                                starttime: Optional[float] = None,
                                timeshift: str = 'first') -> Tuple[int, float]:
        """
        Identifies a cutoff value for the calculation based on the fraction
        of acf samples less than the pure noise std.

        Parameters
        ----------
        delta : int
            The size of the rolling window to use.
        threshold : float or None, optional
            The threshold value to use.  The cutoff is identified based on the
            first acf value where the fraction of samples inside the window
            greater than the std crosses the threshold. Default value is 0.8.
        startindex : int or None, optional
            All acf data starting with this index will be used to compute the
            standard deviation.  This value should be sufficiently large such
            that the acf at that point is pure noise.  Cannot be given with
            starttime.  If neither startindex or starttime are given,
            startindex will be taken as the midpoint of the acf array.
        starttime : float or None, optional
            All acf data starting with this time value will be used to compute
            the standard deviation.  This value should be sufficiently large such
            that the acf at that point is pure noise.  Cannot be given with
            startindex.  If neither startindex or starttime are given,
            startindex will be taken as the midpoint of the acf array.
        timeshift : str, optional
            Indicates where the time value returned is taken from within the
            rolling window: "first" (default), "middle" or "last".
        """

        # Compute the std noise fraction
        f = self.std_noise_fraction(delta, startindex=startindex, starttime=starttime)

        # Identify the first index where f > threshold
        icut: int = np.where(f >= threshold)[0][0]

        # Shift according to the timeshift option
        if timeshift == 'last':
            icut += delta - 1
        elif timeshift == 'middle':
            icut += (delta - 1) // 2
        elif timeshift != 'first':
            raise ValueError('invalid timeshift option')
        
        return icut, self.time[icut]





class GreenKuboKappa(GreenKubo):

    def __init__(self,
                 time: npt.ArrayLike,
                 heatflux: Optional[npt.ArrayLike] = None,
                 acf: Optional[npt.ArrayLike] = None,
                 acfsamples: Union[int, npt.ArrayLike] = 1,
                 temperature: Optional[float] = None,
                 volume: Optional[float] = None):
        """
        Class providing supporting analysis tools for Green-Kubo
        equilibrium calculations of thermal conductivity.

        Parameters
        ----------
        time : array-like
            The values of time for the heatflux or acf values.
        heatflux : array-like or None, optional
            The heat flux values that the auto-correlation function is computed for,
            i.e. <J(t)*J(0)>. Either value or acf must be given.
        acf : array-like or None, optional
            A pre-computed auto-correlation function to evaluate.  Either
            heatflux or acf must be given.
        acfsamples : int or array-like, optional
            If acf is given this indicates how many (relative) samples were
            averaged over for each given acf value.  This is typically either a
            constant (default = 1) if the same number of samples were made for each
            acf, or is a linearly decreasing array if the acf was computed for an
            ongoing known time series value.  This is automatically determined if
            x is given.
        temperature : float or None, optional
            The system temperature to use when computing the thermal conductivity.
            Can be set during init or by setting to the associated class attribute.
        volume : float or None, optional
            The system volume to use when computing the thermal conductivity.
            Can be set during init or by setting to the associated class attribute.
        """
        super().__init__(time, x=heatflux, acf=acf, acfsamples=acfsamples)
        
        self.temperature = temperature
        self.volume = volume

    @classmethod
    def acf_from_lammps(cls,
                        filename: str,
                        name: str,
                        timestep: float,
                        unit: Optional[str] = None,
                        index: int = 2147483647,
                        temperature: Optional[float] = None,
                        volume: Optional[float] = None):
        """
        Read in step and auto-correlation data from LAMMPS fix ave/correlate
        output.

        Parameters
        ----------
        filename : str or Path
            The data file generated by the fix ave/correlate file.
        name : str
            The column name associated with the auto-correlation function data.
        timestep : float
            The timestep used by the LAMMPS simulation, i.e.
            time is TimeDelta * timestep.
        unit : str or None, optional
            The units that the acf data is in. If None (default), no unit
            conversions will be done.
        index : int, optional
            The fix ave/correlate file may have multiple time series runs.
            This allows for the selection of which series to read in.  If
            index is greater than the number of time series runs in the file,
            the final series will be read in (which is what is usually wanted).
            The default value is 2147483647, the 32-bit maximum int, which
            should in a practical sense always be larger than the number of series.
        temperature : float or None, optional
            The system temperature to use when computing the thermal conductivity.
            Can be set during init or by setting to the associated class attribute.
        volume : float or None, optional
            The system volume to use when computing the thermal conductivity.
            Can be set during init or by setting to the associated class attribute.
        """
        obj = super(GreenKuboKappa, cls).acf_from_lammps(
            filename, name, timestep, unit=unit, index=index)
        obj.temperature = temperature
        obj.volume = volume
        
        return obj

    @property
    def heatflux(self) -> Optional[np.ndarray]: 
        """numpy.ndarray or None : The heat flux values"""
        return self.x

    @property
    def temperature(self) -> float: 
        """float : The equilibrium simulation temperature"""
        if self.__temperature is None:
            raise ValueError('temperature not set!')
        return self.__temperature
    
    @temperature.setter
    def temperature(self, val: Optional[float]): 
        if val is not None:
            val = float(val)
            assert val > 0
        self.__temperature = val
    
    @property
    def volume(self) -> float: 
        """float or None: The total system volume"""
        if self.__volume is None:
            raise ValueError('volume not set!')
        return self.__volume

    @volume.setter
    def volume(self, val: Optional[float]): 
        if val is not None:
            val = float(val)
            assert val > 0
        self.__volume = val

    def kappa(self, isJV: bool = False) -> np.ndarray:
        """
        Computes the thermal conductivity (kappa) from the cumulative
        integral of the auto-correlation function. NOTE: temperature
        and volume must be set prior to calling this method!

        Parameters
        ----------
        isJV : bool
            Setting this to True indicates that the heat flux given or used
            to compute given acf values were actually J*V and not J. For 
            instance, the LAMMPS compute heat/flux command directly returns
            J*V values.  This setting simply indicates if the kappa calculation
            multiplies or divides by the volume to get the appropriate value.
        
        Returns
        -------
        kappa : numpy.ndarray
            The thermal conductivity values estimated for all values of the
            cumulative integral of the auto-correlation function.
        """

        # Short names for pretty equations
        T = self.temperature
        V = self.volume
        I = self.integral
        kB = uc.unit['kB']

        # Compute kappa from I, V and T
        if isJV:
            # Divide by volume if flux used was J*V
            kappa = I / (V * T**2 * kB)
        else:
            # Multiply by volume if flux used was J
            kappa = I * V / (T**2 * kB)

        return kappa








class GreenKuboMu(GreenKubo):

    def __init__(self,
                 time: npt.ArrayLike,
                 pressure: Optional[npt.ArrayLike] = None,
                 acf: Optional[npt.ArrayLike] = None,
                 acfsamples: Union[int, npt.ArrayLike] = 1,
                 temperature: Optional[float] = None,
                 volume: Optional[float] = None):
        """
        Class providing supporting analysis tools for Green-Kubo
        equilibrium calculations of viscosity.

        Parameters
        ----------
        time : array-like
            The values of time for the pressure or acf values.
        pressure : array-like or None, optional
            The shear pressure values that the auto-correlation function is
            computed for, i.e. <P(t)*P(0)>. Either value or acf must be given.
        acf : array-like or None, optional
            A pre-computed auto-correlation function to evaluate.  Either
            pressure or acf must be given.
        acfsamples : int or array-like, optional
            If acf is given this indicates how many (relative) samples were
            averaged over for each given acf value.  This is typically either a
            constant (default = 1) if the same number of samples were made for each
            acf, or is a linearly decreasing array if the acf was computed for an
            ongoing known time series value.  This is automatically determined if
            x is given.
        temperature : float or None, optional
            The system temperature to use when computing the thermal conductivity.
            Can be set during init or by setting to the associated class attribute.
        volume : float or None, optional
            The system volume to use when computing the thermal conductivity.
            Can be set during init or by setting to the associated class attribute.
        """
        super().__init__(time, x=pressure, acf=acf, acfsamples=acfsamples)
        
        self.temperature = temperature
        self.volume = volume

    @classmethod
    def acf_from_lammps(cls,
                        filename: str,
                        name: str,
                        timestep: float,
                        unit: Optional[str] = None,
                        index: int = 2147483647,
                        temperature: Optional[float] = None,
                        volume: Optional[float] = None):
        """
        Read in step and auto-correlation data from LAMMPS fix ave/correlate
        output.

        Parameters
        ----------
        filename : str or Path
            The data file generated by the fix ave/correlate file.
        name : str
            The column name associated with the auto-correlation function data.
        timestep : float
            The timestep used by the LAMMPS simulation, i.e.
            time is TimeDelta * timestep.
        unit : str or None, optional
            The units that the acf data is in. If None (default), no unit
            conversions will be done.
        index : int, optional
            The fix ave/correlate file may have multiple time series runs.
            This allows for the selection of which series to read in.  If
            index is greater than the number of time series runs in the file,
            the final series will be read in (which is what is usually wanted).
            The default value is 2147483647, the 32-bit maximum int, which
            should in a practical sense always be larger than the number of series.
        temperature : float or None, optional
            The system temperature to use when computing the thermal conductivity.
            Can be set during init or by setting to the associated class attribute.
        volume : float or None, optional
            The system volume to use when computing the thermal conductivity.
            Can be set during init or by setting to the associated class attribute.
        """
        obj = super(GreenKuboMu, cls).acf_from_lammps(
            filename, name, timestep, unit=unit, index=index)
        obj.temperature = temperature
        obj.volume = volume
        
        return obj

    @property
    def pressure(self) -> Optional[np.ndarray]: 
        """numpy.ndarray : The pressure values"""
        return self.x

    @property
    def temperature(self) -> float: 
        """float : The equilibrium simulation temperature"""
        if self.__temperature is None:
            raise ValueError('temperature not set!')
        return self.__temperature
    
    @temperature.setter
    def temperature(self, val: Optional[float]): 
        if val is not None:
            val = float(val)
            assert val > 0
        self.__temperature = val
    
    @property
    def volume(self) -> float: 
        """float : The total system volume"""
        if self.__volume is None:
            raise ValueError('volume not set!')
        return self.__volume

    @volume.setter
    def volume(self, val: Optional[float]): 
        if val is not None:
            val = float(val)
            assert val > 0
        self.__volume = val

    def mu(self) -> np.ndarray:
        """
        Computes the viscosity (mu) from the cumulative integral of the
        auto-correlation function.  NOTE: temperature and volume must
        be set prior to calling this method!
        
        Returns
        -------
        mu : numpy.ndarray
            The viscosity values estimated for all values of the
            cumulative integral of the auto-correlation function.
        """

        # Short names for pretty equations
        T = self.temperature
        V = self.volume
        I = self.integral
        kB = uc.unit['kB']

        # Calculate mu
        mu = (I * V) / (kB * T)

        return mu
    


class GreenKuboDiffusion(GreenKubo):
        
    def __init__(self,
                 time: npt.ArrayLike,
                 velocity: Optional[npt.ArrayLike] = None,
                 acf: Optional[npt.ArrayLike] = None,
                 acfsamples: Union[int, npt.ArrayLike] = 1):
        """
        Class providing supporting analysis tools for Green-Kubo
        equilibrium calculations of diffusion.

        Parameters
        ----------
        time : array-like
            The values of time for the pressure or acf values.
        velocity : array-like or None, optional
            The atomic velocity values that the auto-correlation function is
            computed for, i.e. <v(t)*v(0)>. Either value or acf must be given.
        acf : array-like or None, optional
            A pre-computed auto-correlation function to evaluate.  Either
            velocity or acf must be given.
        acfsamples : int or array-like, optional
            If acf is given this indicates how many (relative) samples were
            averaged over for each given acf value.  This is typically either a
            constant (default = 1) if the same number of samples were made for each
            acf, or is a linearly decreasing array if the acf was computed for an
            ongoing known time series value.  This is automatically determined if
            x is given.
        """
        super().__init__(time, x=velocity, acf=acf, acfsamples=acfsamples)

    @property
    def velocity(self) -> Optional[np.ndarray]: 
        """numpy.ndarray : The velocity values"""
        return self.x
    
    def D(self):
        """
        Computes the diffusion constant (D) from the cumulative integral of the
        auto-correlation function. 
        
        Returns
        -------
        D : numpy.ndarray
            The diffusion values estimated for all values of the
            cumulative integral of the auto-correlation function.
        """
        # D is simply the integral
        return self.integral