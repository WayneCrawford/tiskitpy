#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Functions to test spectral functions
"""
# from os import system
import unittest
from pathlib import Path
from copy import deepcopy

from obspy.core import Stream, Trace

from tiskitpy.utils import SeisRotate


class TestMethods(unittest.TestCase):
    """
    Test suite
    """
    def setUp(self):
        self.path = Path(__file__).resolve().parent
        self.test_path = self.path / "data" / "utils"

    def test_get_one_trace(self):
        tr = SeisRotate._get_one_trace(
            _quick_stream(['BH1', 'BH2', 'BH3', 'BHZ', 'BHP']), '1')
        self.assertEqual(tr.stats.channel, 'BH1')
        # Stream with missing component should raise IndexError
        with self.assertRaises(IndexError):
            SeisRotate._get_one_trace(_quick_stream(['BH1', 'BH2', 'BH3']), 'Z')
        # Stream with redundant component should raise ValueError
        with self.assertRaises(ValueError):
            SeisRotate._get_one_trace(_quick_stream(['BH1', 'BH2', 'SH2']), '2')

    def test_get_seis_traces(self):
        Z, N, E = SeisRotate._get_seis_traces(
            _quick_stream(['BH1', 'BH2', 'BH3', 'BHZ', 'BHP']))
        self.assertEqual(Z.stats.channel, 'BHZ')
        self.assertEqual(N.stats.channel, 'BH1')
        self.assertEqual(E.stats.channel, 'BH2')

    def test_separate_streams(self):
        seis, other = SeisRotate.separate_streams(
            _quick_stream(['BH1', 'BH2', 'BH3', 'BHZ', 'BHP']))
        self.assertEqual(len(seis), 3)
        self.assertEqual(len(other), 2)
        seis, other = SeisRotate.separate_streams(
            _quick_stream(['BHE', 'BHN', 'BHZ']))
        self.assertEqual(len(seis), 3)
        self.assertIsNone(other)
        # Stream with missing component should raise IndexError
        with self.assertRaises(IndexError):
            SeisRotate.separate_streams(
                _quick_stream(['BH1', 'BH2', 'BH3', 'BHP']))
        # Stream with redundant component should raise ValueError
        with self.assertRaises(ValueError):
            SeisRotate.separate_streams(_quick_stream(['BH1', 'BH2', 'BHZ', 'SH2']))

        
def _quick_stream(chan_list):
    "Stream with given channels"
    stream = Stream()
    for chan in chan_list:
        stream += Trace(header={"channel": chan})
    return stream


def suite():
    return unittest.makeSuite(TestMethods, 'test')


if __name__ == '__main__':
    unittest.main(defaultTest='suite')
