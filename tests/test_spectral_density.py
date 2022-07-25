#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Functions to test spectral functions
"""
import unittest
import inspect
import copy
import logging
from pathlib import Path

import numpy as np

from tiskit import SpectralDensity
from make_test_stream import make_test_stream


class TestMethods(unittest.TestCase):
    """
    Test suite
    """

    def setUp(self):
        self.path = (Path(inspect.getfile(inspect.currentframe()))
                     .resolve().parent)
        self.test_path = self.path / "data" / "spectral_density"
        self.stream, self.sineparms = make_test_stream()
        self.sd = SpectralDensity.from_stream(self.stream, window_s=100.0)

    def test_str(self):
        """Test __str__ function"""
        self.assertEqual(
            self.sd.__str__(),
            "SpectralDensity object:\n"
            "\tchannels=['XX.STA.00.BX1', 'XX.STA.00.BX2', 'XX.STA.00.BX3']\n"
            "\tchannel_units=['Counts', 'Counts', 'Counts']\n"
            "\t8192 frequencies, from 0.0061 to 50Hz\n"
            "\tn_windows=6\n"
            "\twindow_type=prol1pi",
        )

    def test_channels(self):
        """Test channels derived property"""
        self.assertEqual(self.sd.channels, ["XX.STA.00.BX1", "XX.STA.00.BX2",
                                            "XX.STA.00.BX3"])

    def test_freqs(self):
        """Test freqs derived property"""
        stop = self.stream[0].stats.sampling_rate / 2
        num = len(self.sd.freqs)
        self.assertTrue(
            np.all(self.sd.freqs == np.arange(1, num + 1) * stop / num)
        )

    def test_other_properties(self):
        """Test other properties"""
        # window_type
        self.assertEqual(self.sd.window_type, "prol1pi")
        sts = self.sd.starttimes
        self.assertEqual(sts[0], self.stream[0].stats.starttime)
        self.assertEqual(len(sts), 6)
        self.assertEqual(sts[1] - sts[0], 16384 / 100.0)
        self.assertEqual(self.sd.n_windows, 6)

    def test_autospect(self):
        """Test autospect retrieval function"""
        x = self.sd.autospect("XX.STA.00.BX1")
        # Verify that it returns a numpy array of the right size
        self.assertEqual(x.size, 8192)
        self.assertIsInstance(x, np.ndarray)
        # Verify that a bad channel name raises a ValueError
        self.assertRaises(ValueError, self.sd.autospect, *["blah"])

    def test_crossspect(self):
        """Test crossspect retrieval function"""
        x = self.sd.crossspect("XX.STA.00.BX1", "XX.STA.00.BX2")
        # Verify that it returns a numpy array of the right size
        self.assertEqual(x.size, 8192)
        self.assertIsInstance(x, np.ndarray)
        # Verify that the phase of the cross-spect equals 45°
        # (self.chan_amp_phases['BX2'][0].phase)
        sineparm = self.sineparms["XX.STA.00.BX2"][0]
        ifreq = np.argmin(np.abs(self.sd.freqs - sineparm.frequency))
        self.assertAlmostEqual(np.degrees(np.angle(x[ifreq])), sineparm.phase,
                               places=2)
        # Verify that crossspect for same channel is same as autospect
        self.assertTrue(
            np.all(
                self.sd.crossspect("XX.STA.00.BX1", "XX.STA.00.BX1")
                == self.sd.autospect("XX.STA.00.BX1")
            )
        )

    def test_verify_channel(self):
        """Test _verify_channel function"""
        # Verify that a bad channel name raises a ValueError
        self.assertRaises(
            ValueError, self.sd._verify_channel, *["blah", "fake_channel"]
        )
        # Verify that a bad channel type raises a TypeError
        self.assertRaises(
            TypeError, self.sd._verify_channel, *[2, "fake_channel"]
        )
        # Verify that a good channel name returns nothing
        self.assertEqual(
            self.sd._verify_channel("XX.STA.00.BX1", "fake_channel"), None
        )

    def test_replace_channel_name(self):
        """Test replace_channel_name function"""
        x = copy.deepcopy(self.sd)
        x.replace_channel_name("XX.STA.00.BX1", "toto")
        self.assertEqual(x.channels, ["toto", "XX.STA.00.BX2",
                                      "XX.STA.00.BX3"])
        self.assertTrue(
            np.all(x.autospect("toto") == self.sd.autospect("XX.STA.00.BX1"))
        )

    def test_channel_units(self):
        """Test channel_units function"""
        # No response in test case, so units are Counts
        self.assertEqual(self.sd.channel_units("XX.STA.00.BX1"), "Counts")

    def test_units(self):
        """Test units function"""
        # No response in test case, so units are (Counts)^2/Hz
        self.assertEqual(
            self.sd.units("XX.STA.00.BX1", "XX.STA.00.BX2"), "(Counts)^2/Hz"
        )

    def test_coherence(self):
        """Test coherence retrieval function"""
        x = self.sd.coherence("XX.STA.00.BX1", "XX.STA.00.BX2")
        # self.sd.plot_coherences()
        # Verify that coherence at shared frequency is nearly 1.0
        f = 13
        ifreq = np.argmin(np.abs(self.sd.freqs - f))
        self.assertAlmostEqual(x[ifreq], 1.0)

    def test_make_class(self):
        # test window_s > data_len
        self.assertRaises(ValueError, SpectralDensity.from_stream,
                          *[self.stream], **dict(window_s=2000.0))

    def test_spectral_density_values(self):
        delfreq = self.sd.freqs[1] - self.sd.freqs[0]
        # self.sd.plot(overlay=True)
        ihf = 3  # half-width of frequency bins to sum
        # First channel should have peak with amplitude 10 at 13 Hz
        s = self.sd.autospect("XX.STA.00.BX1")
        sp = self.sineparms["XX.STA.00.BX1"][0]
        ifreq = np.argmin(np.abs(self.sd.freqs - sp.frequency))
        self.assertAlmostEqual(
            np.sqrt(np.sum(s[ifreq - ihf: ifreq + ihf]) * delfreq),
            sp.amplitude, delta=0.1)
        # Second channel should have peaks with amplitudes:
        # 2 at 13 Hz, 1 at 10 Hz and 3 at 42 Hz
        s = self.sd.autospect("XX.STA.00.BX2")
        # for f, a in zip((13, 10, 42), (2, 1, 3)):
        for sp in self.sineparms["XX.STA.00.BX2"]:
            ifreq = np.argmin(np.abs(self.sd.freqs - sp.frequency))
            self.assertAlmostEqual(
                np.sqrt(np.sum(s[ifreq - ihf: ifreq + ihf]) * delfreq),
                sp.amplitude, delta=0.1)

    def test_sliding_window(self):
        """Test the _sliding_window() function"""
        # Test a window = data but whos power of two is larger
        logging.disable()
        a = self.stream[0].data  # 1000 s at 100 sps
        logging.disable(logging.NOTSET)
        # Test a window longer than the data
        self.assertRaises(
            ValueError, SpectralDensity._sliding_window, *[a, 200000]
        )
        # Test a window almost as long as the data
        _, nd, _ = SpectralDensity._sliding_window(a, 65536)
        self.assertEqual(nd, 1)
        # Test a window much smaller than the data
        _, nd, _ = SpectralDensity._sliding_window(a, 64)
        self.assertEqual(nd, 1562)
        # Test normal window size
        _, nd, _ = SpectralDensity._sliding_window(a, 10000)
        self.assertEqual(nd, 10)
        _, nd, _ = SpectralDensity._sliding_window(a, 10001)
        self.assertEqual(nd, 9)
        # Test overlap
        _, nd, _ = SpectralDensity._sliding_window(a, 10000, ss=9000)
        self.assertEqual(nd, 11)

    def test_remove_outliers(self):
        """Test the _remove_outliers() method"""
        # Create a 15-element data array with 3 elements whose variance is > 3
        # 2 with variance > 10, 1 with variance > 100
        t = np.arange(0, 100, 0.01)
        amps = [1.05, 1.1, 1.1, 1.1, 1.1, 1.09, 1.11, 1.12, 1.1, 1.13, 1.15, 2.5, 3.5, 20, 200]
        data = np.sin(t) + 1j*np.cos(t)
        A = [am*data for am in amps]
        # for am in amps:
        #     A.append(am*data)
        fts = {'test': np.array(A),
               'ch2': 10*np.array(A)}
        sts = ['test' for x in A]
        self.assertEqual(len(sts), len(amps))
        
        # Test the method
        ft, newsts = SpectralDensity._remove_outliers(fts, sts)
        self.assertEqual(len(newsts), 11)
        ft, newsts = SpectralDensity._remove_outliers(fts, sts, z_threshold=3.5)
        self.assertEqual(len(newsts), 13)
        ft, newsts = SpectralDensity._remove_outliers(fts, sts, z_threshold=2.2)
        self.assertEqual(len(newsts), 9)


def suite():
    return unittest.makeSuite(TestMethods, "test")


if __name__ == "__main__":
    unittest.main(defaultTest="suite")
