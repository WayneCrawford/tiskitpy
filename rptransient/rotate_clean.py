#!env python3
"""
Clean tilt noise from OBS vertical channel using non-deforming rotation
"""
import sys
import math as M
from copy import deepcopy

import numpy as np
import scipy as sp
# import matplotlib as mp
import matplotlib.pyplot as plt
import obspy
from obspy.signal.rotate import rotate2zne
from obspy.core.stream import Stream  # , Trace
from obspy import UTCDateTime
# from obspy.signal import rotate

# from .EQ_template import EQTemplate
from .eq_spans import get_eq_spans

DEBUG = True


def debug(s):
    if DEBUG:
        print(s)


def rotate_clean(stream, excludes=[], horiz_too=False, plot=False,
                 quickTest=False, remove_eq=True, uselogvar=False,
                 verbose=True, filt_band=(0.001, 0.01)):
    """
    Rotates vertical channel to minimize noise

    Arguments:
        stream (Stream): input data, must have a *Z, *[1|N] and *[2|E] channel
        remove_eq (str, True or False): filename of catalog to use to remove
            earthquakes, will download catalog from USGS if True, not remove
            EQs if False
        excludes: list of dictionaries containing time periods to avoid,
                     using "start" and "end" keys
        horiz_too: rotate horizontals also (use True if you think
                      the channels are truly orthogonal)
        plot: Plot comparision of original and rotated vertical
        quickTest: Only run one day's data and do not save results
        uselogvar(bool): use logarithm of variance as metric
        filt_band (tuple): lower, upper frequency limits of band to filter data
            before calculating rotation

    Returns:
        (tuple): 3-tuple containing:
            (Stream): rotated stream
            (float): angle by which Z (or Z-X-Y) was rotated
            (float): azimuth by which Z (or Z-X-Y) was rotated
    """
    if isinstance(remove_eq, str):
        eq_spans = get_eq_spans(stream[0].stats.starttime,
                                  stream[0].stats.endtime, eqfile=remove_eq,
                                  quiet=not verbose)
    if remove_eq is True:
        eq_spans = get_eq_spans(stream[0].stats.starttime,
                                  stream[0].stats.endtime,
                                  quiet=not verbose)
    else:
        eq_spans = None
    # Filter data to tilt noise band for best Angle calc)
    filtstream = stream.copy()
    for interval in excludes:
        filtstream.cutout(UTCDateTime(interval['start']),
                          UTCDateTime(interval['end']))
    filtstream.detrend('demean')
    filtstream.detrend('linear')
    filtstream.filter('lowpass', freq=filt_band[1], corners=4, zerophase=True)
    filtstream.filter('highpass', freq=filt_band[0], corners=4, zerophase=True)

    if verbose:
        print('Creating seisRotate object')
    srData = SeisRotate(filtstream)

    if verbose:
        print('Calculating rotation angle')
    (bestAngle, bestAzimuth) = srData.calc_zrotate_opt(verbose=verbose,
                                                       eq_spans=eq_spans,
                                                       uselogvar=uselogvar)
    if verbose:
        print('Best angle,azimuth is ({:.2f},{:.2f}).  Rotating axes'
              .format(bestAngle, bestAzimuth))

    # Filter data (to make it easier to verify that rotation works
    # on low freq noise)
    viewstream = stream.copy()
    if verbose:
        print('Filtering data for viewing')
    viewstream.filter('lowpass', freq=0.005, corners=4, zerophase=True)
    viewstream.filter('highpass', freq=0.001, corners=4, zerophase=True)

    srData = SeisRotate(stream)
    viewData = SeisRotate(viewstream)

    srData.zrotate(bestAngle, bestAzimuth, horiz_too)
    viewData.zrotate(bestAngle, bestAzimuth)

    if verbose:
        print('Converting back to Stream')
    strm_rot = srData.stream()
    view_rot = viewData.stream()

    # PLOT RESULTS (Z channels)
    trace_view_Z = viewstream.select(component='Z')[0]
    trace_view_rot_Z = view_rot.select(component='Z')[0]
    trace_view_rot_Z.stats.channel = trace_view_rot_Z.stats.channel[0] + \
        'R' + trace_view_rot_Z.stats.channel[2]
    compare_stream = obspy.core.stream.Stream([trace_view_Z,
                                               trace_view_rot_Z])
    if plot:
        compare_stream.plot(equal_scale=True, method='full')

    return strm_rot, bestAngle, bestAzimuth


def rot_tfs(angle, azimuth):
    """
    Return the Z-1 and Z-2 transfer functions equivalent to the given rotation
    
    Designed to be used with the output of rotate_clean()
    I DID THIS ON THE FLY, HAVE NOT VERIFIED THE VALUES

    Arguments:
        angle: rotation angle (degrees)
        azimuth: rotation azimuth (degrees)

    Returns:
        (tuple): 2-tuple containing:
            (float): Z-1 ratio
            (float): Z-2 ration
    """
    
    # Calculate the horizontal to vertical ratio for the given angle
    hratio = np.sin(np.radians(angle))
    Z1_ratio = np.abs(hratio * np.cos(np.radians(azimuth)))
    Z2_ratio = np.abs(hratio * np.sin(np.radians(azimuth)))
    return Z1_ratio, Z2_ratio


class SeisRotate:
    """
    Class to clean tilt noise from OBS vertical channel through
    non-deforming rotation
    """
    def __init__(self, stream, uselogvar=False):
        """
        Create a seisRotate object from a 3-component obsPy Stream
        
        Arguments:
            stream (Stream): 3-component data stream
            uselogvar(bool): use logarithmic variance when searching for
                             best angles

        Channel names must end in Z, N and E or Z, 1, and 2
        Z is up, 1 and 2 are horizontal orthogonal with 2 90° clockwise
        of "1"  (in other words, "1" corresponds to "N" and "2" to "E",
        except that they are not necessarily aligned with geographic
        cardinals)
        """
        self.uselogvar = uselogvar
        # Read in data
        try:
            self.Z = stream.select(component='Z')[0].copy()
        except Exception:
            print('No Z channel found')
            sys.exit(2)
        try:
            self.N = stream.select(component="N")[0].copy()
        except Exception:
            try:
                self.N = stream.select(component="1")[0].copy()
            except Exception:
                print('No Z or "1" channel found')
                sys.exit(2)
        try:
            self.E = stream.select(component="E")[0].copy()
        except Exception:
            try:
                self.E = stream.select(component="2")[0].copy()
            except Exception:
                print('No E or "2" channel found')
                sys.exit(2)

        # Verify all channels have same length
        assert len(self.Z.data) == len(self.N.data), \
            print('Z and N vectors have different lengths')
        assert len(self.Z.data) == len(self.E.data), \
            print('Z and E vectors have different lengths')

        # Verify all channels have same sampling_rate
        self.fs = self.Z.stats.sampling_rate
        assert self.fs == self.E.stats.sampling_rate
        assert self.fs == self.N.stats.sampling_rate

    def __len__(self):
        return len(self.Z.data)

    def plot(self):
        self.stream().plot()

    def copy(self):
        """
        Returns:
            (SeisRotate): Copy of SeisRotate object.
        """
        return deepcopy(self)

    def stream(self):
        """ Return obspy 3-component stream()"""
        return Stream(traces=[self.Z, self.N, self.E])

    def zrotate(self, angle, azimuth, horiz_also=True):
        """
        Rotate a seismometer's vertical component
        
        Arguments:
            angle (float): angle (towards vertical) by which to rotate (degrees)
            azimuth (float): azimuth (CW from N) at which to rotate (degrees)
            horiz_also (bool): rotate horizontal components as well? (to retain
                               orthogonality)

        I would set horiz_also to:
          False for the Guralp CMG-3T, which independently levels the three
                   components
          True for the Trillium T240, which centers without rotation
                   (according to Dr Wieland)
        """
        dip_Z = -90 + angle

        # Calculated rotations for N and E using following logic
        # if azimuth=0 (aligned North), then Ndip=vertang, Edip=0
        # if azimuth=90 (aligned East), then Ndip=0, Edip=vertang
        # if azimuth=180 (aligned South), then Ndip=-vertang, Edip=0
        # if azimuth=270 (aligned West), then Ndip=0, Edip=-vertang
        #
        # Fits Ndip=vertang*cos(azimuth), Edip=vertang*sin(azimuth)
        dip_N = angle * M.cos(np.deg2rad(azimuth))
        dip_E = angle * M.sin(np.deg2rad(azimuth))

        [self.Z.data, N, E] = rotate2zne(
            self.Z.data, azimuth, dip_Z,
            self.N.data, 0, dip_N,
            self.E.data, 90, dip_E)
        if horiz_also:
            self.N.data = N
            self.E.data = E

    def calc_zrotate_opt(self, lowcut=0.001, hicut=0.005, eq_spans=None,
                         uselogvar=None, verbose=False):
        """
        Calculate the Z channel rotation angle that minimizes tilt noise

        Arguments:
            lowcut (float): low passband frequency in which to lood for energy
                             reduction
            hicut (float): high passpand frequency "" "" ""
            eq_spans (TimeSpans): information on time spans to ignore
            uselogvar(bool): use logarithmic variance estimate
            verbose (bool): output information about the angles tested?

        The default (lowcut, hicut) values of (0.001, 0.005) correspond to the
        band where removing tilt noise generally has the greatest effect on
        the Z signal (at the low frequency end of the noise notch)

        You can input the returned angles to rotateSeis.zrotate() to get
        the properly rotated data
        """
        if uselogvar is not None:
            self.uselogvar = uselogvar
        # Filter data, cut off edge effects and detrend
        filt = self.copy()
        filt.Z.filter('bandpass', freqmin=lowcut, freqmax=hicut, corners=5)
        filt.N.filter('bandpass', freqmin=lowcut, freqmax=hicut, corners=5)
        filt.E.filter('bandpass', freqmin=lowcut, freqmax=hicut, corners=5)
        if eq_spans is not None:
            filt.Z= eq_spans.zero(filt.Z)
            filt.N= eq_spans.zero(filt.N)
            filt.E= eq_spans.zero(filt.E)

        # Quick estimate of angles using Z/E and Z/N ratios. DOESN'T HELP: SKIP
        # (startAngle,startAzimuth)=filt._estimateAngles(verbose)
        startAngle, startAzi = (0, 0)

        # Search for the best angles
        angle, azimuth = filt._searchBestAngles(startAngle, startAzi,
                                                verbose)
        return angle, azimuth

    def _estimateAngles(self, verbose=False):
        """
        Estimate how far and in which direction Z is tilted from vertical
        by comparing with N and E

        azimuth is 0 towards N, 90 towards E
        dip is w.r.t. vertical

        Assumes data has already been filtered into the relevant band

        Doesn't seem to work AT ALL! I think there is too much other stuff
        (transients, etc) for this to work.  Also doesn't account for any time
        lags: might work better if ratios calculated from cross-corellation

        Arguments:
            verbose (bool): display extra information?
        """
        if verbose:
            print("Estimating preliminary angle based on signal ratios")
        # s=obspy.core.stream.Stream([self.Z,self.N])
        # s.plot()
        ZoverN = np.divide(self.Z.data, self.N.data)
        ZoverE = np.divide(self.Z.data, self.E.data)
        # Throw out not-a-numbers
        ZoverN = ZoverN[~np.isnan(ZoverN)]
        ZoverE = ZoverE[~np.isnan(ZoverE)]
        # if verbose:
        #   print("    ",np.median(ZoverN),np.median(ZoverE))

        ZNratio = M.copysign(M.sqrt(np.median(ZoverN**2)), np.median(ZoverN))
        ZEratio = M.copysign(M.sqrt(np.median(ZoverE**2)), np.median(ZoverE))
        ZNratio = np.median(ZoverN)
        ZEratio = np.median(ZoverE)
        plt.plot(ZoverN, 'x')

        if verbose:
            print(f'    {ZNratio=}, {ZEratio=}')

        ZNangle = (180/M.pi) * M.asin(ZNratio)
        ZEangle = (180/M.pi) * M.asin(ZEratio)
        azimuth = (180/M.pi) * M.atan(ZNratio / ZEratio)
        azimuth = 90 - (180/M.pi) * M.atan2(ZNratio, ZEratio)
        if verbose:
            print(f'    ZNangle={ZNangle:<10g}, ZEangle={ZEangle:<10g}')
        # The angle from vertical below is just a guess, should developed
        # mathematically (or at least derived empirically from the searched
        # results)
        angle = M.sqrt(ZNangle**2 + ZEangle**2)
        if verbose:
            print('    estAngle= {angle:<10g}, estAzim = {azimuth:<10g}')

        # Sanity check: do the dip equations from zrotate()
        # return ZN & ZE angles??
        dip_N = angle * M.cos(np.deg2rad(azimuth))
        dip_E = angle * M.sin(np.deg2rad(azimuth))
        print(f'{dip_N=}, {dip_E=}')

        return angle, azimuth

    def _searchBestAngles(self, startAngle=0, startAzimuth=0,
                          verbose=False):
        """
        Find best Z rotation angles

        Arguments:
            startAngle (float): starting guess for the angle (degrees)
            startAzimuth (float): starting guess for the azimuth (degrees)
            verbose (bool): display extra information?

        Returns:
            (tuple): 2-tuple containing:
                angle (float): best angle (degrees)
                azimuth (float): best azimuth (degrees)

        Searches for minimum Z energy as function of angle
        """
        if verbose:
            print("Calculating best angle based on variance minimization")
        start_var = self._rotZ_variance([startAngle, startAzimuth])
        if verbose:
            print(f'Starting variance = {start_var}')

        xopt, fopt, iter, funcalls, warnflag, allvecs = sp.optimize.fmin(
                                    func=self._rotZ_variance,
                                    x0=[startAngle, startAzimuth],
                                    disp=verbose,
                                    full_output=True,
                                    retall=True)
        bestAngles = xopt
        if verbose:
            print(f'{xopt=}, {fopt=}, {iter=}, {funcalls=}, {warnflag=}')
            print("    best Angle,Azimuth = {:.2f},{:.2f}"
                  .format(bestAngles[0], bestAngles[1]))
            if self.uselogvar is False:
                print("    variance reduced from {:.2e} to {:.2e} ({:.1f}% lower)"
                    .format(start_var, fopt, 100*(1 - fopt/start_var)))
            if self.uselogvar is False:
                print(" log variance reduced from {:.1f} to {:.1f} ({:.1f}% lower)"
                    .format(start_var, fopt, 100*(1 - 10**(fopt-start_var))))
        if warnflag:
            print(f'{allvecs=}')
            
        return (bestAngles[0], bestAngles[1])

    def _rotZ_variance(self, angles):
        """
        Calculate the variance for a given rotation

        Arguments:
            angles (list): angle, azimuth (in degrees)

        Assumes data are already filtered into relevant band
        """
        A = self.copy()
        A.zrotate(angles[0], angles[1])
        if self.uselogvar is True:
            var = np.log10(np.sum(A.Z.data**2))
        else:
            var = np.sum(A.Z.data**2)
        return var
