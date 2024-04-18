#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Functions to test the lcheapo functions
"""
# from __future__ import (absolute_import, division, print_function,
#                         unicode_literals)
# from future.builtins import *  # NOQA @UnusedWildImport
# 
import unittest
import inspect
from pathlib import Path

from obspy import read

from tiskitpy.fir_corr import fir2caus


class TestFirCorrMethods(unittest.TestCase):
    """
    Test suite for fircorr operations.
    """
    def setUp(self):
        self.path = Path(inspect.getfile(
            inspect.currentframe())).resolve().parent
        self.testing_path = self.path / "data" / "fir_corr"

    def test_fircorr(self):
        """
        Test running fircorr.
        """
        # Read in original and fir-corrected data
        s_orig = read(str(self.testing_path /
                          "example_orig_LSVEI_20150827061217.mseed"), 'MSEED')
        s_fir = read(str(self.testing_path /
                         "example_FIR_LSVEI_20150827061217.mseed"), 'MSEED')

        # Correct orig data using fir_corr routine
        s_fird = fir2caus(s_orig,
                          str(self.testing_path / 'lc2000_fir3_0.json'),
                          2)
        self.assertListEqual(s_fir[0].data.tolist(), s_fird[0].data.tolist())


def suite():
    return unittest.makeSuite(TestFirCorrMethods, 'test')


if __name__ == '__main__':
    unittest.main(defaultTest='suite')
