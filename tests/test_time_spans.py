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

from tiskitpy import TimeSpans


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
        base = UTCDateTime(2011, 11, 1)
        ts = TimeSpans([[base, base+100],
                         [base+1000, base+1100],
                         [base+2000, base + 2100]])
                                
        # Start and End Times before and after time_spans
        self.assertEqual(ts.invert(base-100, base+3000),
                         TimeSpans([[base-100, base],
                                    [base+100, base+1000],
                                    [base+1100, base+2000],
                                    [base+2100, base+3000]]))
        # Start after start of time_spans
        self.assertEqual(ts.invert(base+50, base+3000),
                         TimeSpans([[base+100, base+1000],
                                    [base+1100, base+2000],
                                    [base+2100, base+3000]]))
        self.assertEqual(ts.invert(base+150, base+3000),
                         TimeSpans([[base+150, base+1000],
                                    [base+1100, base+2000],
                                    [base+2100, base+3000]]))
        # End before end of time_spans
        self.assertEqual(ts.invert(base-100, base+2000),
                         TimeSpans([[base-100, base],
                                    [base+100, base+1000],
                                    [base+1100, base+2000]]))
        self.assertEqual(ts.invert(base-100, base+2500),
                          TimeSpans([[base-100, base],
                                    [base+100, base+1000],
                                    [base+1100, base+2000],
                                    [base+2100, base+2500]]))
        # Start after time_spans start and end before time_spans end
        self.assertEqual(ts.invert(base+200, base+2100),
                         TimeSpans([[base+200, base+1000],
                                    [base+1100, base+2000]]))
        # Start and end between two time spans
        self.assertEqual(ts.invert(base+200, base+900),
                         TimeSpans([[base+200, base+900]]))
        # Unspecified Start and End Times: use first end_time as start and
        # last start_time as end
        self.assertEqual(ts.invert(None, None),
                         TimeSpans([[base+100, base+1000],
                                    [base+1100, base+2000]]))
        # Start after time_spans
        self.assertEqual(ts.invert(base+5000, base+6000),
                         TimeSpans([[base+5000, base+6000]]))
        # End before time_spans
        self.assertEqual(ts.invert(base-5000, base-2000),
                         TimeSpans([[base-5000, base-2000]]))
        # ERRORS:
        # Start after end
        with self.assertRaises(ValueError):
            ts.invert(base+1000, base+100)
        # Start equals end
        with self.assertRaises(ValueError):
            ts.invert(base+1000, base+1000)
       
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

    def test_from_eqs(self):
        """Test `from_eqs()` method"""
        eq_spans = TimeSpans.from_eqs(
            UTCDateTime('2013-07-17T00:00'), UTCDateTime('2013-07-28T00:00'),
            minmag=5.5, days_per_magnitude=1.5,
            eq_file = str(self.test_path / 'test_cat.qml'))
        self.assertEqual(eq_spans, TimeSpans(
                [[UTCDateTime('2013-07-15T13:59:04.59'), UTCDateTime('2013-07-17T20:37:43.18')], 
                 [UTCDateTime('2013-07-19T11:40:42.00'), UTCDateTime('2013-07-19T18:52:42.00')],
                 [UTCDateTime('2013-07-20T19:17:10.14'), UTCDateTime('2013-07-23T11:49:42.50')],
                 [UTCDateTime('2013-07-24T03:32:33.59'), UTCDateTime('2013-07-24T17:56:33.59')],
                 [UTCDateTime('2013-07-26T07:07:15.63'), UTCDateTime('2013-07-28T02:20:59.99')]
                ]))

def suite():
    return unittest.makeSuite(TestMethods, 'test')


if __name__ == '__main__':
    unittest.main(defaultTest='suite')
