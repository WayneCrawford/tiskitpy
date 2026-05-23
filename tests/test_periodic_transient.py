#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Functions to test the FIRConverter class
"""
import unittest
import inspect
from pathlib import Path

from obspy.core.stream import Trace
from obspy import UTCDateTime
import numpy as np

from tiskitpy import PeriodicTransient


class TestPTMethods(unittest.TestCase):
    """
    Test suite for PeriodicTransient operations.
    """
    def setUp(self):
        self.path = Path(inspect.getfile(
            inspect.currentframe())).resolve().parent
        # self.testing_path = self.path / "data" / "fir_conv"

    def test_calc_slice_starttime(self, plot=False):
        """
        Test running FIRConverter.
        """
        pt = PeriodicTransient('test', 3600, 1., (-1000, 1000), '2026-10-02')

        trace = Trace(np.ones(40), {'starttime': UTCDateTime('2026-10-02')})        
        self.assertEqual(pt._transient_offset(trace), 0.)
        self.assertEqual(pt._calc_slice_starttime(trace),
                         trace.stats.starttime)

        trace = Trace(np.ones(40),
                      {'starttime': UTCDateTime('2026-10-01T23:30:00')})
        self.assertEqual(pt._transient_offset(trace), 1800)
        self.assertEqual(pt._calc_slice_starttime(trace),
                         UTCDateTime('2026-10-01T23:40:00'))

def suite():
    return unittest.makeSuite(TestPTMethods, 'test')


if __name__ == '__main__':
    unittest.main(defaultTest='suite')
