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

from tiskitpy import SpectralDensity, TimeSpans
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
        stream2 = self.stream.copy()
        trace = stream2.select(channel='*3')[0]
        trace.stats['clean_sequence'] = ['ROT']
        self.sd_tagged = SpectralDensity.from_stream(stream2, window_s=100.)
        

    def test_str(self):
        """Test __str__ function"""
        self.assertEqual(
            self.sd.__str__(),
            "SpectralDensity object:\n"
            "\tids=['XX.STA.00.BX1', 'XX.STA.00.BX2', 'XX.STA.00.BX3', 'XX.STA.00.BDH']\n"
            "\tchannel_units=['Counts', 'Counts', 'Counts', 'Counts']\n"
            "\t8192 frequencies, from 0.0061 to 50Hz\n"
            "\tn_windows=6\n"
            "\twindow_type=prol1pi",
        )

    def test_ids(self):
        """Test derived property "ids" """
        self.assertEqual(self.sd.ids, ["XX.STA.00.BX1", "XX.STA.00.BX2",
                                       "XX.STA.00.BX3",  "XX.STA.00.BDH"])
        self.assertEqual(self.sd_tagged.ids, ["XX.STA.00.BX1", "XX.STA.00.BX2",
                                              "XX.STA.00-ROT.BX3",  "XX.STA.00.BDH"])
        self.assertEqual(self.sd_tagged.seed_ids, ["XX.STA.00.BX1", "XX.STA.00.BX2",
                                                   "XX.STA.00.BX3",  "XX.STA.00.BDH"])
        self.assertEqual(self.sd_tagged.seed_id("XX.STA.00-ROT.BX3"),  "XX.STA.00.BX3")
        self.assertEqual(self.sd_tagged.seed_id("XX.STA.00.BX2"),  "XX.STA.00.BX2")

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
        sineparm = self.sineparms["XX.STA.00.BX2"][0]
        ifreq = np.argmin(np.abs(self.sd.freqs - sineparm.frequency))
        self.assertAlmostEqual(np.degrees(np.angle(x[ifreq])), sineparm.phase,
                               places=2)
        # Verify that crossspect for same channel is same as autospect
        chan = "XX.STA.00.BX1"
        self.assertTrue(np.all(self.sd.crossspect(chan, chan)
                        == self.sd.autospect(chan)))

    def test_channel_id(self):
        """Test channel_id() function"""
        # Verify that a bad channel id raises a ValueError
        with self.assertRaises(ValueError):
            self.sd.channel_id("blah")
        # Verify that a bad channel type raises a TypeError
        with self.assertRaises(TypeError):
            self.sd.channel_id(2)
        # Verify that a good channel name returns nothing
        self.assertEqual(self.sd.channel_id("XX.STA.00.BX1"), "XX.STA.00.BX1")
        self.assertEqual(self.sd.channel_id("*BX1"), "XX.STA.00.BX1")
        # Verify that multiple fitting channel name raises a ValueError
        with self.assertRaises(ValueError):
            self.sd.channel_id("*.BX*")

    def test_replace_channel_id(self):
        """Test replace_channel_id() method"""
        x = copy.deepcopy(self.sd)
        x.replace_channel_id("XX.STA.00.BX1", "toto")
        self.assertEqual(x.ids, ["toto", "XX.STA.00.BX2",
                                 "XX.STA.00.BX3", 'XX.STA.00.BDH'])
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
        logging.disable()
        npts = 100000
        logging.disable(logging.NOTSET)
        # Window longer than the data
        self.assertEqual(SpectralDensity._sliding_window(npts, 200000), [])
        # Window almost as long as the data
        offsets = SpectralDensity._sliding_window(npts, 65536)
        self.assertEqual(len(offsets), 1)
        # Window much smaller than the data
        offsets = SpectralDensity._sliding_window(npts, 64)
        self.assertEqual(len(offsets), 1562)
        # Normal window size
        offsets = SpectralDensity._sliding_window(npts, 10000)
        self.assertEqual(len(offsets), 10)
        offsets = SpectralDensity._sliding_window(npts, 10001)
        self.assertEqual(len(offsets), 9)
        # Overlap
        offsets = SpectralDensity._sliding_window(npts, 10000, ss=9000)
        self.assertEqual(len(offsets), 11)

    def test_make_windows(self):
        """Test the _make_windows() function"""
        # Test a window = data but whos power of two is larger
        logging.disable()
        tr = self.stream[0]  # 1000 s at 100 sps
        logging.disable(logging.NOTSET)
        cls = SpectralDensity
        tr_start = tr.stats.starttime
        tr_end = tr.stats.endtime
        sr = 100
        wind_s = 100
        ws = sr*wind_s  # Standard window length in samples
        ws_toolong = 200000  # Window length > data

        # Test error cases
        # no TimeSpan within the data
        ts = TimeSpans(start_times=[tr_end+wind_s], end_times=[tr_end+10*wind_s])
        a, sts = cls._make_windows(tr, ws, ws, "hanning", None, ts)
        # starttimes outside of trace times
        starttimes = [tr_start + x for x in (-1, 100, 200)]
        with self.assertRaises(ValueError):
            cls._make_windows(tr, ws, ws, "hanning", starttimes, None)
        starttimes = [tr_start + x for x in (0, 100, 2000)]
        with self.assertRaises(ValueError):
            cls._make_windows(tr, ws, ws, "hanning", starttimes, None)
        
        # Test non-errors
        starttimes = [tr_start + x for x in (0, 100, 200, 300)]
        ts = TimeSpans(start_times=[tr_start], end_times=[tr_end+1])
        a, sts = cls._make_windows(tr, ws, ws, "hanning", None, None)
        self.assertEqual(a.shape, (10, ws))
        self.assertEqual(a.shape[0], len(sts))
        a, sts = cls._make_windows(tr, ws, ws, "hanning", None, ts)
        self.assertEqual(a.shape, (10, ws))
        self.assertEqual(a.shape[0], len(sts))
        a, sts = cls._make_windows(tr, ws, ws, "hanning", starttimes, None)
        self.assertEqual(a.shape, (len(starttimes), ws))
        self.assertEqual(a.shape[0], len(sts))
        # time_spans outside of specified range (just skips bad times)
        ts = TimeSpans(start_times=[tr_start+wind_s], end_times=[tr_end+wind_s])
        a, sts = cls._make_windows(tr, ws, ws, "hanning", None, ts)
        self.assertEqual(a.shape, (9, ws))
        ts = TimeSpans(start_times=[tr_start-wind_s], end_times=[tr_end-wind_s+1])
        a, sts = cls._make_windows(tr, ws, ws, "hanning", None, ts)
        self.assertEqual(a.shape, (9, ws))
        # Window longer than the data
        self.assertEqual(cls._make_windows(tr, ws_toolong, ws_toolong,
                         "hanning",None, None), (None, None))
       
        # Test other errors
        # Bad window_taper name
        with self.assertRaises(ValueError):
            cls._make_windows(tr, ws, ws, "harpo", None, None)
        # Both starttimes and time_spans specified
        with self.assertRaises(RuntimeError):
            cls._make_windows(tr, ws, ws, "hanning", starttimes, ts)
       

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
        self.assertEqual(len(newsts), 13)
        ft, newsts = SpectralDensity._remove_outliers(fts, sts, z_threshold=3.5)
        self.assertEqual(len(newsts), 15)
        ft, newsts = SpectralDensity._remove_outliers(fts, sts, z_threshold=2.2)
        self.assertEqual(len(newsts), 9)


def suite():
    return unittest.makeSuite(TestMethods, "test")


if __name__ == "__main__":
    unittest.main(defaultTest="suite")
