#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Functions to test the lcheapo functions
"""
# from os import system
import unittest
import inspect
from pathlib import Path

from obspy.core import UTCDateTime
from obspy.core.stream import read as stream_read

from tiskit import TimeSpans


class TestMethods(unittest.TestCase):
    """
    Test suite
    """
    def setUp(self):
        self.path = Path(inspect.getfile(
            inspect.currentframe())).resolve().parent
        self.test_path = self.path / "data" / "time_spans"

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


def suite():
    return unittest.makeSuite(TestMethods, 'test')


if __name__ == '__main__':
    unittest.main(defaultTest='suite')
