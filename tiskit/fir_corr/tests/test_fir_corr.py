#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Functions to test the lcheapo functions
"""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
from future.builtins import *  # NOQA @UnusedWildImport

import os
import unittest
import inspect
import sys
from obspy import read

# For now, just point to fircorr module
sys.path.append('..')
from fir2caus import fir2caus


class TestFirCorrMethods(unittest.TestCase):
    """
    Test suite for fircorr operations.
    """
    def setUp(self):
        self.path = os.path.dirname(os.path.abspath(inspect.getfile(
            inspect.currentframe())))
        self.testing_path = os.path.join(self.path, "data")

    def test_fircorr(self):
        """
        Test running fircorr.
        """
        # Read in original and fir-corrected data
        s_orig = read(os.path.join(self.testing_path,
                                   "example_orig_LSVEI_20150827061217.mseed"),
                      'MSEED')
        s_fir = read(os.path.join(self.testing_path,
                                  "example_FIR_LSVEI_20150827061217.mseed"),
                     'MSEED')

        # Correct orig data using fir_corr routine
        s_fird = fir2caus(s_orig,
                          os.path.join(self.testing_path,
                                       'lc2000_fir3_0.json'),
                          2)
        self.assertListEqual(s_fir[0].data.tolist(), s_fird[0].data.tolist())


def suite():
    return unittest.makeSuite(TestFirCorrMethods, 'test')


if __name__ == '__main__':
    unittest.main(defaultTest='suite')
