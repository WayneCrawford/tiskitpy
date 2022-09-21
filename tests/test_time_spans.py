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

    def test_overlaps(self):
        """
        Test overlaps in a TimeSpan object
        """
        start_times = [UTCDateTime(2011, 11, 1),
                       UTCDateTime(2011, 11, 2),
                       UTCDateTime(2011, 12, 1)]
        end_times = [x + 100000 for x in start_times]
        ts = TimeSpans(start_times=start_times, end_times=end_times)
        self.assertEqual(len(ts), 2)
        self.assertEqual(
            ts, TimeSpans(start_times=[UTCDateTime(2011, 11, 1),
                                       UTCDateTime(2011, 12, 1)],
                          end_times=[UTCDateTime(2011, 11, 2) + 100000,
                                     UTCDateTime(2011, 12, 1) + 100000]))

    def test_creator(self):
        """
        Test TimeSpan creator
        """
        start_times = [UTCDateTime(2011, 11, 1), UTCDateTime(2011, 12, 1)]
        end_times = [x + 86400 for x in start_times]
        spans = [[x, x+86400] for x in start_times]
        ts1 = TimeSpans(start_times=start_times, end_times=end_times)
        ts2 = TimeSpans(spans)
        ts3 = TimeSpans(start_times=['2011-11-01', '2011-12-01'],
                        end_times=['2011-11-02', '2011-12-02'])
        ts4 = TimeSpans([['2011-11-01', '2011-11-02'],
                         ['2011-12-01', '2011-12-02']])
        self.assertEqual(ts1, ts2)
        self.assertEqual(ts2, ts3)
        self.assertEqual(ts3, ts4)
        # Does this belong here?
        self.assertEqual(ts1.start_times, start_times)
        self.assertEqual(ts1.end_times, end_times)
        self.assertEqual(ts1.spans, spans)

        # CHECK ERRORS
        # Input value that can't be converted to a UTCDateTime
        with self.assertRaises(TypeError):
            TimeSpans([['2011-11-01', '2011-11-02'],
                       ['2011-12-01', 'BLAH']])
        # Non-increasing span
        with self.assertRaises(ValueError):
            TimeSpans([['2011-11-01', '2011-11-02'],
                       ['2011-12-04', '2011-12-02']])
            TimeSpans(starttimes=['2011-11-01', '2011-12-04'],
                      end_times=['2011-11-02', '2011-12-02'])
        # Different-length start_times and end_times
        with self.assertRaises(ValueError):
            TimeSpans(start_times=['2011-11-01'],
                      end_times=['2011-11-02', '2011-12-02'])

    def test_properties(self):
        """
        Test TimeSpan object, including combining overlaps
        """
        start_times = [UTCDateTime(2011, 11, 1),
                       UTCDateTime(2011, 11, 2),
                       UTCDateTime(2011, 12, 1)]
        end_times = [x + 1000 for x in start_times]
        ts = TimeSpans(start_times=start_times, end_times=end_times)
        self.assertEqual(start_times, ts.start_times)
        self.assertEqual(end_times, ts.end_times)
        self.assertEqual(ts.spans,
                         [[x, y] for x, y in zip(start_times, end_times)])

    def test_invert(self):
        """
        Test TimeSpan object's invert() method'.
        """
        start_times = [UTCDateTime(2011, 11, 1),
                       UTCDateTime(2011, 11, 2),
                       UTCDateTime(2011, 12, 1)]
        end_times = [x + 1200 for x in start_times]  # 20-minute windows
        time_spans = TimeSpans(start_times=start_times, end_times=end_times)
        ts_start = UTCDateTime(2011, 10, 30)
        ts_end = UTCDateTime(2011, 12, 2)
        
        # Should raise error when time series start time is AFTER the first
        # time span
        with self.assertRaises(ValueError):
            time_spans.invert(start_times[1], ts_end)

        # Should raise error when given time series end time is BEFORE the end
        # of the last time span
        with self.assertRaises(ValueError):
            time_spans.invert(ts_start, start_times[-1])
        
        # now check that invert() does what is expected
        new_starts= [ts_start]
        new_starts.extend(time_spans.end_times)
        new_ends = time_spans.start_times
        new_ends.append(ts_end)
        self.assertEqual(time_spans.invert(ts_start, ts_end),
                         TimeSpans(start_times=new_starts, end_times=new_ends))

    def test_time_spans_zero_interp(self):
        """
        Test TimeSpan object's zero() and interp() methods.
        """
        stream = stream_read(str(self.test_path / 'XS.S10D.LH.mseed'))
        st = UTCDateTime(2016, 12, 5, 6)
        et = UTCDateTime(2016, 12, 5, 12)
        ts = TimeSpans([[st, et]])  # , save_eq_file=False, quiet=True)
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
