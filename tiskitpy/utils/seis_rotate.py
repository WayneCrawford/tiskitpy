#!env python3
"""
Clean tilt noise from OBS vertical channel using non-deforming rotation

Is there any reason why rotate_clean isn't just rotate_calc + rotate_apply?
"""
import math as M
from copy import deepcopy

import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
from obspy.signal.rotate import rotate2zne
from obspy.core.stream import Stream  # , Trace

from ..logger import init_logger
from .stream_synchronize import stream_synchronize

logger = init_logger()


class SeisRotate:
    """
    Class to clean tilt noise from OBS vertical channel through
    non-deforming rotation
    """

    def __init__(self, stream, uselogvar=False, max_reject_sync=0.01):
        """
        Create a seisRotate object from a 3-component obsPy Stream

        Arguments:
            stream (Stream): 3+ component data stream
            uselogvar(bool): use logarithmic variance when searching for
                             best angles
            max_reject_sync(): max_reject value for stream_synchronize()

        Channel names must end in Z, N and E or Z, 1, and 2
        Z is up, 1 and 2 are horizontal orthogonal with 2 90° clockwise
        of "1"  (in other words, "1" corresponds to "N" and "2" to "E",
        except that they are not necessarily aligned with geographic
        cardinals)
        """
        stream = stream_synchronize(stream, max_reject_sync)
        self.uselogvar = uselogvar
        self.Z, self.N, self.E = SeisRotate._get_seis_traces(stream)
        self.fs = self.Z.stats.sampling_rate
        # Verify that all channels have same length
        if not len(self.Z.data) == len(self.N.data):
            raise ValueError("Z and N vectors have different lengths")
        if not len(self.Z.data) == len(self.E.data):
            raise ValueError("Z and E vectors have different lengths")
        # Verify that all channels have same sampling_rate
        if not self.Z.stats.sampling_rate == self.E.stats.sampling_rate:
            raise ValueError("Z and E sampling rates are different")
        if not self.Z.stats.sampling_rate == self.N.stats.sampling_rate:
            raise ValueError("Z and N sampling rates are different")

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
        """Return obspy 3-component stream()"""
        return Stream(traces=[self.Z, self.N, self.E])

    def zrotate(self, angle, azimuth, horiz_also=True):
        """
        Rotate a seismometer's vertical component

        Arguments:
            angle (float): angle towards vertical by which to rotate (degrees)
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
            self.Z.data,
            azimuth,
            dip_Z,
            self.N.data,
            0,
            dip_N,
            self.E.data,
            90,
            dip_E,
        )
        if horiz_also:
            self.N.data = N
            self.E.data = E

    def calc_zrotate_opt(
        self,
        lowcut=0.001,
        hicut=0.005,
        ignore_spans=None,
        uselogvar=None
    ):
        """
        Calculate the Z channel rotation angle that minimizes tilt noise

        Arguments:
            lowcut (float): low passband frequency in which to evaluate
                            variance reduction
            hicut (float): high passpand frequency "" "" ""
            ignore_spans (:class:`TimeSpans``): time spans to ignore
            uselogvar(bool): use logarithmic variance estimate

        Returns: (tuple)
            angle (float): tilt angle (degrees, 0 means no tilt correction)
            azimuth (float): tilt azimuth (degrees)
            var_red (float): obtained reduction in variance (0 to 1)

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
        filt.Z.filter("bandpass", freqmin=lowcut, freqmax=hicut, corners=5)
        filt.N.filter("bandpass", freqmin=lowcut, freqmax=hicut, corners=5)
        filt.E.filter("bandpass", freqmin=lowcut, freqmax=hicut, corners=5)
        if ignore_spans is not None:
            filt.Z = ignore_spans.zero(filt.Z)
            filt.N = ignore_spans.zero(filt.N)
            filt.E = ignore_spans.zero(filt.E)

        # Quick estimate of angles using Z/E and Z/N ratios. DOESN'T HELP: SKIP
        startAngle, startAzi = (0, 0)

        # Search for the best angles
        angle, azimuth, var_red = filt._searchBestAngles(startAngle, startAzi)
        return angle, azimuth, var_red

    def _estimateAngles(self):
        """
        Estimate how far and in which direction Z is tilted from vertical
        by comparing with N and E

        azimuth is 0 towards N, 90 towards E
        dip is w.r.t. vertical

        Assumes data has already been filtered into the relevant band

        Doesn't seem to work AT ALL! I think there is too much other stuff
        (transients, etc) for this to work.  Also doesn't account for any time
        lags: might work better if ratios calculated from cross-corellation

        """
        logger.debug("Estimating preliminary angle based on signal ratios")
        ZoverN = np.divide(self.Z.data, self.N.data)
        ZoverE = np.divide(self.Z.data, self.E.data)
        # Throw out not-a-numbers
        ZoverN = ZoverN[~np.isnan(ZoverN)]
        ZoverE = ZoverE[~np.isnan(ZoverE)]

        ZNratio = M.copysign(M.sqrt(np.median(ZoverN**2)), np.median(ZoverN))
        ZEratio = M.copysign(M.sqrt(np.median(ZoverE**2)), np.median(ZoverE))
        ZNratio = np.median(ZoverN)
        ZEratio = np.median(ZoverE)
        plt.plot(ZoverN, "x")

        logger.debug(f"    {ZNratio=}, {ZEratio=}")

        ZNangle = (180 / M.pi) * M.asin(ZNratio)
        ZEangle = (180 / M.pi) * M.asin(ZEratio)
        azimuth = (180 / M.pi) * M.atan(ZNratio / ZEratio)
        azimuth = 90 - (180 / M.pi) * M.atan2(ZNratio, ZEratio)
        logger.debug(f"    ZNangle={ZNangle:<10g}, ZEangle={ZEangle:<10g}")
        # The angle from vertical below is just a guess, should developed
        # mathematically (or at least derived empirically from the searched
        # results)
        angle = M.sqrt(ZNangle**2 + ZEangle**2)
        logger.debug("    estAngle= {angle:<10g}, estAzim = {azimuth:<10g}")

        # Sanity check: do the dip equations from zrotate()
        # return ZN & ZE angles??
        dip_N = angle * M.cos(np.deg2rad(azimuth))
        dip_E = angle * M.sin(np.deg2rad(azimuth))
        logger.info(f"{dip_N=}, {dip_E=}")

        return angle, azimuth

    def _searchBestAngles(self, startAngle=0, startAzimuth=0):
        """
        Find best Z rotation angles

        Arguments:
            startAngle (float): starting guess for the angle (degrees)
            startAzimuth (float): starting guess for the azimuth (degrees)

        Returns:
            (tuple): 2-tuple containing:
                angle (float): best angle (degrees, 0 to 90)
                azimuth (float): best azimuth (degrees, 0 to 360)
                var_red (float): variance reduction (0 to 1)

        Searches for minimum Z energy as function of angle
        """
        logger.debug("Calculating best angle based on variance minimization")
        start_var = self._rotZ_variance([startAngle, startAzimuth])
        logger.debug(f"Starting variance = {start_var}")

        xopt, fopt, iter, funcalls, warnflag, allvecs = sp.optimize.fmin(
            func=self._rotZ_variance,
            x0=[startAngle, startAzimuth],
            disp=False, full_output=True, retall=True)
        bestAngle, bestAzimuth = xopt
        if bestAngle < 0:
            bestAngle = -bestAngle
            bestAzimuth += 180
        bestAzimuth %= 360.
        if bestAngle > 90:
            logger.warning(f'bestAngle > 90! ({bestAngle})')

        logger.debug(f"{xopt=}, {fopt=}, {iter=}, {funcalls=}, {warnflag=}")

        if self.uselogvar is False:
            var_red = 1 - fopt / start_var
            logger.debug("    variance reduced from {:.2e} to {:.2e} ({:.1f}% lower)"
                        .format(start_var, fopt, 100 * var_red))
        else:
            var_red = 1 - 10 ** (fopt - start_var)
            logger.debug("    log variance reduced from {:.1f} to {:.1f} ({:.1f}% lower)"
                        .format(start_var, fopt, 100 * var_red))
        if warnflag:
            logger.debug(f"optimization warning flag: {allvecs=}")

        return (bestAngle, bestAzimuth, var_red)

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

    @staticmethod
    def separate_streams(stream):
        """
        Separates input stream into seismometer and "other"

        There must be one channel name ending in Z, one ending in N or 1 and
        one ending in E or 2

        Arguments:
            stream (Stream): 3+ component data stream
        Returns:
            seis_stream (Stream): stream containt 3C seismometer channels
            other_stream (Stream): all other channels, or None
        """
        Z, N, E = SeisRotate._get_seis_traces(stream)
        if len(stream) == 3:
            return stream.copy(), None
        other_stream = stream.copy()
        seis_stream = Stream(traces=(Z, N, E))
        other_stream.remove(Z)
        other_stream.remove(N)
        other_stream.remove(E)
        return seis_stream, other_stream

    @staticmethod
    def _get_seis_traces(stream):
        Z = SeisRotate._get_one_trace(stream, "Z")
        try:
            N = SeisRotate._get_one_trace(stream, "N")
        except IndexError:
            N = SeisRotate._get_one_trace(stream, "1")
        try:
            E = SeisRotate._get_one_trace(stream, "E")
        except IndexError:
            E = SeisRotate._get_one_trace(stream, "2")
        return Z, N, E

    @staticmethod
    def _get_one_trace(stream, component):
        sel_stream = stream.select(component=component)
        if len(sel_stream) == 0:
            raise IndexError(f'no "{component}" channel found')
        elif len(sel_stream) > 1:
            raise ValueError('more than one "{}" channel: {}'.format(
                component, [x.id for x in sel_stream]))
        return sel_stream[0].copy()
