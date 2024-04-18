#!env python3
""" Decimate data to a lower rate"""

###############################################################################
# Note: resample() takes too long, even on 1.25 sps data
#       Could (should?) write:
#         -  resample extension that uses FIR filter (min phase or zero phase,
#            use scipy.signal.firwin to calculate FIR)
#         - decimate extension that uses FIR filter (min or zero phase,
#           scipy.signal.firwin)
#       Both should return the FIR filter, which could be stuffed into the
#       instruments' response information
#       scipy doesn't seem to have a minimum phase FIR method: could use
#       Scherbaum method to convert zero phase to minimum phase

# import obspy
# from obspy.clients.filesystem import sds
import time
# from copy import deepcopy
from dataclasses import dataclass
from inspect import getfile, currentframe
from pathlib import Path
import warnings
import fnmatch
# import logging

from numpy import prod
from obspy.core.stream import Stream, Trace
from obspy.core.trace import Stats
from obspy.core.inventory import (FIRResponseStage,
                                  CoefficientsTypeResponseStage)

from .fir_filter import FIRFilter
from ..logger import init_logger

logger = init_logger()


@dataclass
class Decimator:
    """
    Class to decimate obspy stream intelligenty

    Can also update the station inventory with the used filter responses

    Args:
        decimates (list): list of decimation factorst to use (integers between
                        2 and 7, will be applied in order)
        verbose (bool): Be chatty
    """
    decimates: list
    verbose: bool = False

    # def __post_init__(self):
    #     if self.verbose is True:
    #         logger.setLevel(logging.INFO)
    #     else:
    #         logger.setLevel(logging.WARN)

    @property
    def decimation_factor(self):
        """Total decimation (product of `decimates`)"""
        return prod(self.decimates)

    def decimate(self, data, keep_dtype=True):
        """
        Apply decimator to data

        Args:
            data (:class:`obspy.Stream` or :class:`obspy.Trace`): waveform data
            keep_dformat (bool): force output data format to be the same as
                the input data_format
        """
        self.keep_dtype = keep_dtype
        if isinstance(data, Stream):
            return self._run_stream(data)
        if isinstance(data, Trace):
            return self._run_trace(data)
        else:
            raise TypeError("data is not a Stream or Trace")

    def update_inventory(
        self, inv, st=None, normalize_firs=False, quiet=False,
        inv_force_overwrite=False
    ):
        """
        Update inventory for channels found in stream

        Args:
            inv (:class:`obspy.core.events.inventory.Inventory`): inventory
            st (:class:`obspy.core.stream.Stream`): data stream (used to
                determine net/sta/chan/loc codes).  If None, will update
                every channel
            normalize_firs (bool): normalize FIR channels that aren't already
            inv_force_overwrite (bool): If the target channel seed_id already
                exists (and it's not the source channel seed_id), overwrite it

        Returns:
            obspy inventory with new channels
        """
        inv = inv.copy()  # Don't overwrite input inventory
        if st is not None:
            stats_list = [tr.stats for tr in st]
        else:
            stats_list = [
                Stats(
                    header=dict(
                        network=net.code,
                        station=sta.code,
                        location=ch.location_code,
                        channel=ch.code,
                        sampling_rate=ch.sample_rate,
                    )
                )
                for net in inv
                for sta in net
                for ch in sta
            ]
        for stats in stats_list:
            new_cha = self._duplicate_channel(
                inv,
                network=stats.network,
                station=stats.station,
                location=stats.location,
                channel=stats.channel,
                normalize_firs=normalize_firs,
                insert_in_inv=False 
            )
            if not new_cha.sample_rate == stats.sampling_rate:
                raise ValueError(
                    "data and metadata sampling rates are different!")
            self._modify_chan(
                new_cha, net=stats.network, sta=stats.station, quiet=quiet
            )
            seed_id = f'{stats.network}.{stats.station}.{stats.location}{new_cha.code}'
            if new_cha.code == stats.channel:
                # Do not overwrite mother channel
                logger.error('New channel has same code as source channel '
                             f'({seed_id}), rejected!')
                continue
            elif len(inv.select(network=stats.network, station=stats.station,
                                location=stats.location, channel=new_cha.code)) > 0:
                boilerplate = 'New channel has same seed_id as an existing channel'
                if inv_force_overwrite is True:
                    logger.warning(f'{boilerplate} ({seed_id}), overwriting!')
                    inv = inv.remove(network=stats.network, station=stats.station,
                                     location=stats.location, channel=new_cha.code)
                else:
                    logger.warning(f'{boilerplate} ({seed_id}), rejected!')
                    continue
            # Append new channel to inv
            appended = False
            for net in inv:
                if not net.code == stats.network:
                    continue
                for sta in net:
                    if sta.code == stats.station:
                        sta.channels.append(new_cha)
                        appended = True
                        break
            if appended is not True:
                raise ValueError(f'Did not find a home for {seed_id}')        
        return inv

    def update_inventory_from_nslc(self, inv, network="*", station="*",
                                   channel="*", location="*", quiet=False,
                                   normalize_firs=False):
        """
        Update inventory based on network, station, channel and location codes

        Args:
            inv (:class:`obspy.core.events.inventory.Inventory`): inventory
            network (str): FDSN network code
            station (str): station code
            channel (str): channel code
            location (str): FDSN location code
            quiet (bool): Do not say what channel(s) was changed
            normalize_firs (bool): normalize FIR channesl that aren't already

        Returns:
            newinv: inv (:class:`obspy.core.events.inventory.Inventory`):
                updated inventory

        network, station, channel and location codes are '*' by default,
        may include wildcards as specified in obspy.core.stream.Stream.select
        """
        inv = inv.copy()  # Don't overwrite input inventory
        inv_selected = inv.select(
            network=network,
            station=station,
            location=location,
            channel=channel,
        )
        for net in inv_selected:
            for sta in net:
                old_sta = sta.copy()
                for cha in old_sta:
                    new_cha = self._duplicate_channel(
                        inv,
                        net.code,
                        sta.code,
                        cha.location_code,
                        cha.code,
                        normalize_firs=normalize_firs,
                    )
                    self._modify_chan(
                        new_cha, net=network, sta=station, quiet=quiet
                    )
        return inv

    @staticmethod
    def get_band_code(in_band_code, sample_rate):
        """
        Return the channel band code based on an input band_code

        :param in_band_code: input band code ('B' if cutoff > 10s,
            'S' otherwise)
        :param sample_rate: input sample rate (sps)
        """
        if len(in_band_code) != 1:
            raise ValueError(f'"{in_band_code=}" is more than one character')
        if in_band_code in "FCHBMLVURPTQ":
            if sample_rate >= 1000:
                return "F"
            elif sample_rate >= 250:
                return "C"
            elif sample_rate >= 80:
                return "H"
            elif sample_rate >= 10:
                return "B"
            elif sample_rate > 1:
                return "M"
            elif sample_rate > 0.3:
                return "L"
            elif sample_rate > 0.03:
                return "V"
            elif sample_rate > 0.003:
                return "U"
            elif sample_rate >= 0.0001:
                return "R"
            elif sample_rate >= 0.00001:
                return "P"
            elif sample_rate >= 0.000001:
                return "T"
            else:
                return "Q"
        elif in_band_code in "GDES":
            if sample_rate >= 1000:
                return "G"
            elif sample_rate >= 250:
                return "D"
            elif sample_rate >= 80:
                return "E"
            elif sample_rate >= 10:
                return "S"
            else:
                warnings.warn("Short period sensor sample rate < 10 sps")
                return "X"
        else:
            raise TypeError(f'Unknown band base code: "{in_band_code}"')

    def _get_new_band_code(self, in_band_code, in_sample_rate,
                           out_sample_rate):
        """
        Return the channel band code

        :param in_band_code: input band code
        :param in_sample_rate: input sample rate (sps)
        :param out_sample_rate: output sample rate (sps)
        """
        if in_band_code != self.get_band_code(in_band_code, in_sample_rate):
            warnings.warn(
                f"Input band code ({in_band_code}) does not match "
                f" input sampling rate ({in_sample_rate})"
            )
        return self.get_band_code(in_band_code, out_sample_rate)

    def _modify_chan(self, cha, net="", sta="", normalize_firs=False,
                     quiet=False):
        """
        modify reponse and name of a channel to correspond to decimation

        Args:
            cha (:class:`obspy.core.inventory.channel`): original channel
            net (str): network code (just for clearer progress printout)
            sta (str): station code (just for clearer progress printout)
            normalize_firs (bool): normalizes any FIR channel that isn't
                already
        """
        old_seed_id = ".".join([net, sta, cha.location_code, cha.code])
        # old_cha = cha.copy()
        input_sample_rate = cha.sample_rate
        self._add_instrument_response(cha, input_sample_rate)
        self._change_chan_loc(cha, input_sample_rate)
        cha.sample_rate /= self.decimation_factor
        seed_id = ".".join([net, sta, cha.location_code, cha.code])
        logger.info("channel modified from "
                    f"{old_seed_id} ({input_sample_rate:g} sps) "
                    f"to {seed_id} ({cha.sample_rate:g} sps) ")

    @staticmethod
    def _normalize_firs(cha):
        """
        Verifies &, if needed, normalizes channel's FIR or Coefficients coeffs

        Args:
            cha (:class:`obspy.core.inventory.channel`): channel
        """
        for stg in cha.response.response_stages:
            if isinstance(stg, FIRResponseStage):
                if stg.symmetry == "NONE":
                    coeff_sum = sum(stg.coefficients)
                elif stg.symmetry == "EVEN":
                    coeff_sum = 2 * sum(stg.coefficients)
                elif stg.symmetry == "ODD":
                    coeff_sum = (
                        2 * sum(stg.coefficients[:-1]) + stg.coefficients[-1]
                    )
                else:
                    raise ValueError(
                        f"Unknown FIR coefficient symmetry: {stg.symmetry}"
                    )
                if abs(coeff_sum - 1) > 0.01:
                    logger.info(f"DECIMATOR: Sum of FIR coeffs = "
                                f"{coeff_sum}, normalizing")
                    stg.coefficients = [x / coeff_sum
                                        for x in stg.coefficients]
            elif isinstance(stg, CoefficientsTypeResponseStage):
                if sum(stg.denominator) == 0:
                    coeff_sum = sum(stg.numerator)
                    if coeff_sum == 0:
                        continue
                else:
                    coeff_sum = sum(stg.numerator) / sum(stg.denominator)
                # descriptor = "Coeff numerator"
                if abs(coeff_sum - 1) > 0.01:
                    logger.info(
                        f"DECIMATOR: sum(numerator coeffs)/sum(denom coeffs)"
                        f" = {coeff_sum}, normalizing")
                    stg.numerator = [x / coeff_sum for x in stg.numerator]

    def _change_chan_loc(self, cha, in_sr, avoid_codes=[]):
        """
        Change channel (or location) code in place, according to decimation

        Arguments:
            cha (Inventory.Channel) channel to be modified (in place)
            in_sr (float): input sample rate
            avoid_codes (list): channel:location codes to avoid
        """
        cha.code, cha.location_code = self._get_chan_loc(
            cha.code, cha.location_code, in_sr, avoid_codes
        )

    def _get_chan_loc(self, cha_code, loc_code, in_sr, avoid_codes=[]):
        """
        Get new channel (or loc) codes from old codes and decimation factor

        Arguments:
            cha_code: input channel code
            loc_code: input location code
            in_sr: input sampling rate
            avoid_codes (list): channel:location codes to avoid

        Returns:
            tuple:
                new_chan_code (str)
                new_loc_code (str)
        """
        if self.decimation_factor == 1:
            raise ValueError("Decimation == 1!")
        out_sr = in_sr / self.decimation_factor
        new_band_code = self._get_new_band_code(cha_code[0], in_sr, out_sr)
        new_cha_code = new_band_code + cha_code[1:]
        if new_cha_code == cha_code:
            # Have to change location code
            try:
                new_loc_code = f"{int(loc_code) + 1:02d}"
            except Exception:
                new_loc_code = "01"
        else:
            new_loc_code = loc_code
        # If chan_loc exists, increment loc_code until we find an empty one
        while True:
            for c in avoid_codes:
                if (
                    new_cha_code == c.split(":")[0]
                    and new_loc_code == c.split(":")[1]
                ):
                    new_loc_code = f"{int(new_loc_code) + 1:02d}"
                    continue
            break
        return new_cha_code, new_loc_code

    def _add_instrument_response(self, cha, input_sample_rate):
        """
        Append  decimation object's instrument response to an existing
        channel's response
        """
        stage_number = cha.response.response_stages[-1].stage_sequence_number
        for d in self.decimates:
            fir_filter = FIRFilter.from_SAC(d)
            stage_number += 1
            cha.response.response_stages.append(
                fir_filter.to_obspy(
                    input_sample_rate, stage_number,
                    cha.response.response_stages[-1].output_units)
            )
            input_sample_rate /= d
        try:
            cha.response.recalculate_overall_sensitivity()
        except Exception:
            i_u = cha.response.instrument_sensitivity.input_units
            if i_u.upper() == 'PA':
                cha.response.instrument_sensitivity.input_units = "M/S"
                cha.response.recalculate_overall_sensitivity()
                cha.response.instrument_sensitivity.input_units = i_u

    def _run_stream(self, stream):
        """
        Decimate obspy Stream
        """
        if not isinstance(stream, Stream):
            raise ValueError("input stream is not an obspy Stream!")
        st = stream.copy()
        newtr = []
        for tr in st.traces:
            newtr.append(self._run_trace(tr))
        st.traces = newtr
        logger.info("New data has {} samples".format([tr.data.size
                                                      for tr in st]))
        st.verify()
        return st

    def _run_trace(self, trace):
        """
        Decimate obspy Trace
        """
        if not isinstance(trace, Trace):
            raise ValueError("input trace is not an obspy Trace!")
        dtype = trace.data.dtype
        tr = trace.copy()
        sr = tr.stats.sampling_rate
        logger.info("Decimating data from {:g} to {:g} Hz ({:d}x)... "
                    .format(sr, sr / self.decimation_factor,
                            self.decimation_factor))
        tic = time.time()
        for d in self.decimates:
            fir_filter = FIRFilter.from_SAC(d)
            tr.data = fir_filter.convolve(tr.data)
            tr.decimate(d, no_filter=True)
        if self.keep_dtype is True and not tr.data.dtype == dtype:
            tr.data = tr.data.astype(dtype)
        tr.stats.channel, tr.stats.location = self._get_chan_loc(
            tr.stats.channel,
            tr.stats.location,
            tr.stats.sampling_rate * self.decimation_factor)
        logger.info("Took {:.1f} seconds".format(time.time() - tic))
        return tr

    def _read_FIR_SAC(decimation):
        """
        Returns a SAC FIR filter
        """
        base_dir = Path(getfile(currentframe())).resolve().parent
        fbase = f"dec{decimation:d}"
        filename = base_dir.glob(fbase)
        if not len(filename) == 1:
            raise NameError('SAC filter file "{fbase}" not found')
        with open(filename) as f:
            f.readline()
            A = f.readline().split()
            divisor = float(A[0])
            nCoeffs = int(A[1])
            coeffs = []
            for line in f.readlines():
                coeffs = coeffs + [float(x) / divisor for x in line.split()]
        if not len(coeffs) == nCoeffs:
            raise RuntimeError(
                "FIR length is different from that stated: "
                "{} != {}".format(len(coeffs), nCoeffs)
            )
        # SAC only saves right half of symmetrical odd filter
        delay = len(coeffs) - 1
        coeffs = coeffs[-1:0:-1] + coeffs
        return FIRFilter(coeffs, delay, 1, "")

    # The following is a combo of obspy's inventory, Network and Station select
    # methods, but which returns the chosen station and channel, not a copy
    def _duplicate_channel(
        self, inv, network, station, location, channel, normalize_firs=False,
        insert_in_inv=True
    ):
        """
        Return Station & Channel matching network, station, location, channel

        Arguments:
            network (str): network code
            station (str): station code
            location (str): location code
            channel (str): channel code
            normalize_firs (bool): normalizes any FIR channel that isn't
                already

        Returns:
            (Channel): duplicate of found channel, ready to modify
        """
        stations = []
        channels = []
        for net in inv.networks:
            # skip if any given criterion is not matched
            if not fnmatch.fnmatch(net.code.upper(), network.upper()):
                continue
            for sta in net.stations:
                if not fnmatch.fnmatch(sta.code.upper(), station.upper()):
                    continue
                for cha in sta.channels:
                    if not fnmatch.fnmatch(
                        cha.location_code.upper(), location.upper()
                    ):
                        continue
                    if not fnmatch.fnmatch(cha.code.upper(), channel.upper()):
                        continue
                    stations.append(sta)
                    channels.append(cha)
        if len(stations) == 0:
            raise ValueError(f'trace waveformid "{network}.{station}.'
                             f'{location}.{channel}" not found in inventory')
            return None, None
        if len(stations) > 1:
            raise ValueError(
                "{:d} station-channels matched {}.{}.{}.{}".format(
                    len(stations), network, station, location, channel
                )
            )
        if normalize_firs is True:
            self._normalize_firs(channels[0])
        new_cha = channels[0].copy()
        if insert_in_inv is True:
            stations[0].channels.append(new_cha)
        return new_cha
