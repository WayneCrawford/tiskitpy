#!env python3
"""
Clean tilt noise from OBS vertical channel using non-deforming rotation

Is there any reason why rotate_clean isn't just rotate_calc + rotate_apply?
"""
import logging

import numpy as np
from obspy.core.stream import Stream
from obspy import UTCDateTime

from .time_spans import TimeSpans
from .utils import SeisRotate


class CleanRotator:
    """
    Clean tilt noise from OBS vertical channel using non-deforming rotation

    Because earthquakes can swamp the noise-based variance, downloads a list
    of earthquakes from the time period and only uses windows outside of the
    earthquakes’ influence (using .TimeSpans.from_eqs).
    Saves the earthquake file locally to speed up future runs

    Args:
        stream (Stream): input data, must have a \\*Z, \\*[1|N] and \\*[2|E]
            channel
        avoid_spans (:class:`TimeSpans`): timespans to avoid
        plot (bool): Plot comparision of original and rotated vertical
        quickTest (bool): Only run one day's data and do not save results
        remove_eq (str, True or False): filename of catalog to use to remove
                earthquakes, will download catalog from USGS if True, not
                remove EQs if False
        uselogvar(bool): use logarithm of variance as metric
        filt_band (tuple): lower, upper frequency limits of band to filter data
                before calculating rotation
        save_eq_file (bool): Passed onto TimeSpans.from_eqs()
    Attributes:
        angle (float): angle by which Z (or Z-X-Y) was rotated
        azimuth (float): azimuth by which Z (or Z-X-Y) was rotated
    """

    def __init__(self, stream, avoid_spans=None, plot=False, quickTest=False,
                 remove_eq=True, uselogvar=False, verbose=True,
                 filt_band=(0.001, 0.01), save_eq_file=True):
        """
        Calculate rotation angles needed to minimize noise on vertical channel

        """
        ignore_spans = self._make_eq_spans(
            remove_eq, stream[0].stats, verbose, save_eq_file
        )
        if avoid_spans is not None:
            ignore_spans += avoid_spans
        filtstream = self._filtstream(stream, filt_band)
        srData = SeisRotate(filtstream)
        (ang, azi) = srData.calc_zrotate_opt(
            verbose=verbose, ignore_spans=ignore_spans, uselogvar=uselogvar
        )
        if verbose:
            logging.info(f"Best angle, azimuth is ({ang:.2f}, {azi:.2f})")
        self.angle = ang
        self.azimuth = azi
        if plot:
            self._plot_filtered_stream(stream, filt_band)

    def __str__(self):
        return "CleanRotator: angle, azimuth = {:.2f}, {:.1f} degrees".format(
            self.angle, self.azimuth
        )

    def _filtstream(self, stream, filt_band):
        """Filter data to tilt noise band for best Angle calc"""
        filtstream = stream.copy()
        filtstream.detrend("demean")
        filtstream.detrend("linear")
        filtstream.filter(
            "lowpass", freq=filt_band[1], corners=4, zerophase=True
        )
        filtstream.filter(
            "highpass", freq=filt_band[0], corners=4, zerophase=True
        )
        return filtstream

    def _make_eq_spans(self, remove_eq, stats, verbose, save_eq_file):
        if isinstance(remove_eq, str):
            return TimeSpans.from_eqs(
                stats.starttime, stats.endtime, quiet=not verbose,
                eq_file=remove_eq, save_eq_file=save_eq_file)
        elif remove_eq is True:
            return TimeSpans.from_eqs(
                stats.starttime, stats.endtime, quiet=not verbose,
                save_eq_file=save_eq_file)
        return None

    def _plot_filtered_stream(self, stream, filt_band):
        viewstream = stream.copy()
        viewstream.filter(
            "lowpass", freq=filt_band[1], corners=4, zerophase=True
        )
        viewstream.filter(
            "highpass", freq=filt_band[0], corners=4, zerophase=True
        )
        viewData = SeisRotate(viewstream)
        viewData.zrotate(self.angle, self.azi)
        view_rot = viewData.stream()
        # PLOT RESULTS (Z channels)
        trace_view_Z = viewstream.select(component="Z")[0]
        trace_view_rot_Z = view_rot.select(component="Z")[0]
        rotZchan = trace_view_rot_Z.stats.channel
        rotZchan = rotZchan[0] + "R" + rotZchan[2]
        compare_stream = Stream([trace_view_Z, trace_view_rot_Z])
        compare_stream.plot(equal_scale=True, method="full")

    def apply(self, stream, horiz_too=False):
        """
        Rotates vertical channel to minimize noise

        Arguments:
            stream (Stream): data, must have \\*Z, \\*[1|N] and \\*[2|E]
                channels
            horiz_too: (bool) rotate horizontals also (use if you believe
                channels are truly orthogonal, probably a bad idea anyway
                as long as we use a 2-value rotation)
        Returns:
            strm_rot (Stream): rotated stream
        """
        seis_stream, other_stream = SeisRotate.separate_streams(stream)
        srData = SeisRotate(stream)
        srData.zrotate(self.angle, self.azimuth, horiz_too)
        if other_stream is None:
            return srData.stream()
        else:
            return srData.stream() + other_stream

    def tfs(self):
        """
        Return the Z-1 and Z-2 transfer functions equal to the given rotation

        Designed to be used with the output of rotate_clean()
        I DID THIS ON THE FLY, HAVE NOT VERIFIED THE VALUES

        Returns:
            (tuple): 2-tuple containing:
                (float): Z-1 ratio
                (float): Z-2 ration
        """
        # Calculate the horizontal to vertical ratio for the given angle
        hratio = np.sin(np.radians(self.angle))
        Z1_ratio = np.abs(hratio * np.cos(np.radians(self.azimuth)))
        Z2_ratio = np.abs(hratio * np.sin(np.radians(self.azimuth)))
        return Z1_ratio, Z2_ratio


def rotate_clean(stream, avoid_spans=None, horiz_too=False, plot=False,
                 quickTest=False, remove_eq=True, uselogvar=False,
                 verbose=True, filt_band=(0.001, 0.01)):
    """
    Rotates vertical channel to minimize noise

    See CleanRotator.__init__() for arguments

    Returns:
        (tuple): 3-tuple containing:
            (Stream): rotated stream
            (float): angle by which Z (or Z-X-Y) was rotated
            (float): azimuth by which Z (or Z-X-Y) was rotated
    """
    obj = CleanRotator(stream, avoid_spans, plot, quickTest, remove_eq,
                       uselogvar, verbose, filt_band)
    return obj.apply(stream), obj.rot_angle, obj.rot_azimuth


# def extract_corr_z(evstream, tf_name):
#     """
#     Return a corrected stream from ATACR EventStream
#
#     Args:
#         evstream (:class:`obstools.atacr.EventStream`):
#         tf_name (str): transfer function key
#     Returns:
#         outstream (Stream): corrected Z stream
#     """
#     stream = evstream.sth.select(component='Z').copy()
#     stream[0].data = evstream.correct[tf_name].flatten()
#     return stream
