#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Functions to test spectral functions
"""
# from os import system
import unittest
from pathlib import Path

from obspy.core import Stream, Trace

# sys.path.append("..")

from tiskit.utils.seis_rotate import SeisRotate as SR


class TestMethods(unittest.TestCase):
    """
    Test suite
    """
    def setUp(self):
        self.path = Path(__file__).resolve().parent
        self.test_path = self.path / "data" / "utils"

    def _quick_stream(self, chan_list):
        "Stream with given channels"
        stream = Stream()
        for chan in chan_list:
            stream += Trace(header={"channel": chan})
        return stream

    def test_get_one_trace(self):
        tr = SR._get_one_trace(
            self._quick_stream(['BH1', 'BH2', 'BH3', 'BHZ', 'BHP']), '1')
        self.assertEqual(tr.stats.channel, 'BH1')
        # Stream with missing component should raise IndexError
        self.assertRaises(IndexError, SR._get_one_trace,
                          *[self._quick_stream(['BH1', 'BH2', 'BH3']), 'Z'])
        # Stream with redundant component should raise ValueError
        self.assertRaises(ValueError, SR._get_one_trace,
                          *[self._quick_stream(['BH1', 'BH2', 'SH2']), '2'])

    def test_get_seis_traces(self):
        Z, N, E = SR._get_seis_traces(
            self._quick_stream(['BH1', 'BH2', 'BH3', 'BHZ', 'BHP']))
        self.assertEqual(Z.stats.channel, 'BHZ')
        self.assertEqual(N.stats.channel, 'BH1')
        self.assertEqual(E.stats.channel, 'BH2')

    def test_separate_streams(self):
        seis, other = SR.separate_streams(self._quick_stream(['BH1', 'BH2',
                                                              'BH3', 'BHZ',
                                                              'BHP']))
        self.assertEqual(len(seis), 3)
        self.assertEqual(len(other), 2)
        seis, other = SR.separate_streams(self._quick_stream(['BHE', 'BHN',
                                                              'BHZ']))
        self.assertEqual(len(seis), 3)
        self.assertIsNone(other)
        # Stream with missing component should raise IndexError
        self.assertRaises(IndexError, SR.separate_streams,
                          *[self._quick_stream(['BH1', 'BH2', 'BH3', 'BHP'])])
        # Stream with redundant component should raise ValueError
        self.assertRaises(ValueError, SR.separate_streams,
                          *[self._quick_stream(['BH1', 'BH2', 'BHZ', 'SH2'])])


def suite():
    return unittest.makeSuite(TestMethods, 'test')


if __name__ == '__main__':
    unittest.main(defaultTest='suite')
