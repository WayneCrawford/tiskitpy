#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Functions to test the lcheapo functions
"""
# from os import system
import unittest
import filecmp
import inspect
import difflib
from pathlib import Path
import pickle
import datetime

from obspy.core import UTCDateTime
from obspy import read_inventory
from obspy.core.stream import read as stream_read
from matplotlib import pyplot as plt
import numpy as np

from rptransient import TimeSpans, rotate_clean

class TestMethods(unittest.TestCase):
    """
    Test suite
    """
    def setUp(self):
        self.path = Path(inspect.getfile(
            inspect.currentframe())).resolve().parent
        self.test_path = self.path / "data"

    def assertTextFilesEqual(self, first, second, msg=None):
        with open(first) as f:
            str_a = f.read()
        with open(second) as f:
            str_b = f.read()

        if str_a != str_b:
            first_lines = str_a.splitlines(True)
            second_lines = str_b.splitlines(True)
            delta = difflib.unified_diff(
                first_lines, second_lines,
                fromfile=first, tofile=second)
            message = ''.join(delta)

            if msg:
                message += " : " + msg

            self.fail("Multi-line strings are unequal:\n" + message)

    def assertBinFilesEqual(self, first, second, msg=None):
        """ Compares two binary files """
        self.assertTrue(filecmp.cmp(first, second))

    def test_time_spans(self):
        """
        Test TimeSpan object.
        """
        start_times = [UTCDateTime(2011, 11, 1),
                       UTCDateTime(2011, 11, 2),
                       UTCDateTime(2011, 12, 1)]
        end_times = [x + 100000 for x in start_times]
        ts = TimeSpans(start_times, end_times)
        self.assertEqual(len(ts), 2)
        self.assertEqual(ts,
                         TimeSpans([UTCDateTime(2011, 11, 1),
                                    UTCDateTime(2011, 12, 1)],
                                   [UTCDateTime(2011, 11, 2) + 100000,
                                    UTCDateTime(2011, 12, 1) + 100000]))

    def test_time_spans_zero_interp(self):
        """
        Test TimeSpan object's zero() method'.
        """
        stream = stream_read(str(self.test_path / 'XS.S10D.LH.mseed'))
        st = UTCDateTime(2016, 12, 5, 6)
        et = UTCDateTime(2016, 12, 5, 12)
        ts = TimeSpans([st], [et])
        zeroed = ts.zero(stream)
        self.assertEqual(zeroed[0].trim(st, et).data[0], 0)
        interped = ts.interp(stream)
        # (stream+interped).select(component='Z').plot()
        self.assertAlmostEqual(
            interped.select(component='Z')[0].trim(st, et).data[0], 213.70,
            delta=0.01)

    def test_rotate_clean(self):
        """
        Test rotate_clean function.
        """
        # this file only has 3 minutes and an EQ, is it enough?
        stream = stream_read(str(self.test_path / 'XS.S10D.LH.mseed'))
        rot_stream, bestAngle, bestAzimuth = rotate_clean(stream, verbose=False) 
        self.assertAlmostEqual(bestAngle, -0.18, delta=0.01)       
        self.assertAlmostEqual(bestAzimuth, 61.67, delta=0.01)       
        

def suite():
    return unittest.makeSuite(TestMethods, 'test')


if __name__ == '__main__':
    unittest.main(defaultTest='suite')
