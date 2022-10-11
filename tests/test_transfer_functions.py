#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Functions to test spectral functions
"""
import unittest
# from pathlib import Path

import numpy as np
from obspy.core import UTCDateTime
from obspy.core.stream import Stream, Trace

from tiskit import SpectralDensity, TransferFunctions
from make_test_stream import make_test_stream


class TestMethods(unittest.TestCase):
    """
    Test suite
    """

    def setUp(self):
        # self.path = (Path(inspect.getfile(inspect.currentframe()))
        #              .resolve().parent)
        # self.test_path = self.path / "data" / "transfer_functions"
        self.stream, self.sineparms = make_test_stream()
        sd = SpectralDensity.from_stream(self.stream, window_s=100.0)
        self.xf = TransferFunctions(sd, "XX.STA.00.BX1")

    def make_stream(
        self,
        freqs,
        chan_amp_phases,
        net="XX",
        sta="STA",
        data_len=1000.0,
        sampling_rate=100,
        loc="00",
        starttime=UTCDateTime(2020, 1, 1),
    ):
        stream = Stream()
        header = dict(
            sampling_rate=sampling_rate,
            network=net,
            station=sta,
            location=loc,
            starttime=starttime,
        )
        for chan in chan_amp_phases.keys():
            header["channel"] = chan
            amps = chan_amp_phases[chan][0]
            phases = chan_amp_phases[chan][1]
            t = np.arange(0, data_len, 1.0 / sampling_rate)
            data = amps[0] * np.sin(
                t * 2 * np.pi * freqs[0] + np.radians(phases[0])
            )
            for f, a, p in zip(freqs[1:], amps[1:], phases[1:]):
                data += a * np.sin(t * 2 * np.pi * f + np.radians(p))
            stream += Trace(data, header)
        return stream

    def test_str(self):
        """Test __str__ function"""
        self.assertEqual(
            self.xf.__str__(),
            "TransferFunctions object:\n"
            "\tinput_channel='XX.STA.00.BX1'\n"
            "\toutput_channels=['XX.STA.00.BX2', 'XX.STA.00.BX3', 'XX.STA.00.BDH']\n"
            "\tnoise_channels=['output', 'output', 'output']\n"
            "\tn_windows=6",
        )

    def test_freqs(self):
        """Test freqs derived property"""
        stop = self.stream[0].stats.sampling_rate / 2
        num = len(self.xf.freqs)
        self.assertTrue(
            np.all(self.xf.freqs == np.arange(1, num + 1) * stop / num)
        )

    def test_coh_signif(self):
        """Test coherence significance levels"""
        self.assertAlmostEqual(self.xf.coh_signif(), 0.88113, places=4)
        self.assertAlmostEqual(self.xf.coh_signif(0.5), 0.541196, places=4)

    def test_channels(self):
        """Test input_channel and output_channel derived properties"""
        self.assertEqual(self.xf.input_channel, "XX.STA.00.BX1")
        self.assertEqual(
            self.xf.output_channels, ["XX.STA.00.BX2", "XX.STA.00.BX3", 'XX.STA.00.BDH']
        )

    def test_other_properties(self):
        """Test other properties"""
        # window_type
        self.assertEqual(self.xf.input_units, "Counts")
        self.assertEqual(self.xf.n_windows, 6)
        self.assertEqual(self.xf.noise_channels, ["output", "output", "output"])
        self.assertEqual(self.xf.noise_channel("XX.STA.00.BX2"), "output")
        # Not much of a test
        self.assertEqual(self.xf.output_units("XX.STA.00.BX2"), "Counts")

    def test_frfs(self):
        """Test frequency response function retrieval and returned values"""
        x = self.xf.frf("XX.STA.00.BX2")
        # Verify that it returns a numpy array of the right size
        self.assertEqual(x.size, 8192)
        self.assertIsInstance(x, np.ndarray)
        f, a_in = 13, .33*1  # Input channel (BX1) Frequency, amplitude
        ifreq = np.argmin(np.abs(self.xf.freqs - f))
        # output channel BX2
        a_out, p_out = .33*.9, 0  # output channel (BX2) amplitude, phase
        pr = np.radians(p_out)
        self.assertAlmostEqual(
            x[ifreq],
            (a_out / a_in) * (np.cos(pr) + 1j * np.sin(pr)),
            delta=0.1,
        )
        # output channel BX3
        x = self.xf.frf("XX.STA.00.BX3")
        a_out, p_out = .33*1e-4, 0  # output channelamplitude, phase
        pr = np.radians(p_out)
        self.assertAlmostEqual(
            x[ifreq],
            (a_out / a_in) * (np.cos(pr) - 1j * np.sin(pr)),
            delta=0.1,
        )

        # Verify that a bad channel name raises a ValueError
        self.assertRaises(ValueError, self.xf.frf, *["blah"])
        # Because there is no transfer function, values() and
        # values_wrt_counts() should be the same
        self.assertTrue(
            np.all(
                self.xf.frf("XX.STA.00.BX2")
                == self.xf.frf_wrt_counts("XX.STA.00.BX2")
            )
        )


def suite():
    return unittest.makeSuite(TestMethods, "test")


if __name__ == "__main__":
    unittest.main(defaultTest="suite")
