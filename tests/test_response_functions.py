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

from tiskitpy import SpectralDensity, ResponseFunctions
from make_test_stream import make_test_stream


class TestMethods(unittest.TestCase):
    """
    Test suite
    """

    def setUp(self):
        self.stream, self.sineparms = make_test_stream()
        sd = SpectralDensity.from_stream(self.stream, window_s=100.0)
        self.xf = ResponseFunctions(sd, "XX.STA.00.BX1")

    def test_str(self):
        """Test __str__ function"""
        self.assertEqual(
            self.xf.__str__(),
            "ResponseFunctions object:\n"
            "  input_channel_id='XX.STA.00.BX1'\n"
            "  output_channel_ids=['XX.STA.00.BX2', 'XX.STA.00.BX3', 'XX.STA.00.BDH']\n"
            "  noise_channels=['output', 'output', 'output']\n"
            "  n_windows=6",
        )

    def test_properties(self):
        """Test derived properties"""
        self.assertEqual(self.xf.input_channel_id, "XX.STA.00.BX1")
        self.assertEqual(self.xf.output_channel_ids,
                         ["XX.STA.00.BX2", "XX.STA.00.BX3", 'XX.STA.00.BDH'])
        self.assertEqual(self.xf.input_clean_sequence, [])
        self.assertEqual(self.xf.input_units, "Counts")
        self.assertEqual(self.xf.n_windows, 6)
        self.assertEqual(self.xf.noise_channels, ["output", "output", "output"])
        self.assertEqual(self.xf.noise_channel("XX.STA.00.BX2"), "output")
        # Not much of a test
        self.assertEqual(self.xf.output_units("XX.STA.00.BX2"), "Counts")
        stop = self.stream[0].stats.sampling_rate / 2
        num = len(self.xf.freqs)
        self.assertTrue(
            np.all(self.xf.freqs == np.arange(1, num + 1) * stop / num)
        )
    def test_coh_signif(self):
        """Test coherence significance levels"""
        self.assertAlmostEqual(self.xf.coh_signif(), 0.88113, places=4)
        self.assertAlmostEqual(self.xf.coh_signif(0.5), 0.541196, places=4)

    def test_values(self):
        """Test frequency response function retrieval and returned values"""
        x = self.xf.value("XX.STA.00.BX2")
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
        x = self.xf.value("XX.STA.00.BX3")
        a_out, p_out = .33*1e-4, 0  # output channelamplitude, phase
        pr = np.radians(p_out)
        self.assertAlmostEqual(
            x[ifreq],
            (a_out / a_in) * (np.cos(pr) - 1j * np.sin(pr)),
            delta=0.1,
        )

        # Verify that a bad channel name raises a ValueError
        self.assertRaises(ValueError, self.xf.value, *["blah"])
        # Because there is no transfer function, values() and
        # values_wrt_counts() should be the same
        self.assertTrue(
            np.all(
                self.xf.value("XX.STA.00.BX2")
                == self.xf.value_wrt_counts("XX.STA.00.BX2")
            )
        )


def suite():
    return unittest.makeSuite(TestMethods, "test")


if __name__ == "__main__":
    unittest.main(defaultTest="suite")
