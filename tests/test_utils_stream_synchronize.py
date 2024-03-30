#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Functions to test spectral functions
"""
import unittest
from pathlib import Path
import inspect

from obspy.core.stream import Stream

from tiskitpy.utils import stream_synchronize
from make_test_stream import make_test_stream

class TestMethods(unittest.TestCase):
    """
    Test suite
    """

    def setUp(self):
        self.stream, _ = make_test_stream()

    def test_stream_synchronize(self):
        """Test stream_synchronize()"""
        # Verify that a good stream is returned as before
        srate = self.stream[0].stats.sampling_rate
        out_stream = stream_synchronize(self.stream)
        self.assertEqual(out_stream, self.stream)

        # Verify that a desynchronized stream is returned corrected
        test_stream = self.stream.copy()
        n_cut = 10
        test_stream[0].data = test_stream[0].data[n_cut:]
        test_stream[0].stats.starttime += n_cut/srate
        out_stream = stream_synchronize(test_stream)
        self.assertNotEqual(test_stream, out_stream)
        for in_tr, out_tr in zip(self.stream, out_stream):
            self.assertEqual(len(in_tr), len(out_tr)+10)
            self.assertAlmostEqual(
                out_tr.stats.starttime - in_tr.stats.starttime,
                10/in_tr.stats.sampling_rate,
                delta=0.5/srate)

        # Verify that a stream that is too desynchronized returns None
        test_stream = self.stream.copy()
        datalen = len(test_stream[0]) 
        test_stream[0].data = test_stream[0].data[int(datalen*.6):]
        self.assertIsNone(stream_synchronize(test_stream))

        # Verify that a non-overlapping stream exits gracefully (returns None)
        test_stream = self.stream.copy()
        one_third = int(len(test_stream[0])/3)
        test_stream[0].data = test_stream[0].data[:one_third]
        test_stream[1].data = test_stream[1].data[-one_third:]
        self.assertIsNone(stream_synchronize(test_stream))


def suite():
    return unittest.makeSuite(TestMethods, "test")


if __name__ == "__main__":
    unittest.main(defaultTest="suite")
