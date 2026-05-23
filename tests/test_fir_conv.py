#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Functions to test the FIRConverter class
"""
import unittest
import inspect
from pathlib import Path

from obspy import read
from matplotlib import pyplot as plt
import numpy as np

from tiskitpy import FIRConverter


class TestFirConvMethods(unittest.TestCase):
    """
    Test suite for FIRConverter operations.
    """
    def setUp(self):
        self.path = Path(inspect.getfile(
            inspect.currentframe())).resolve().parent
        self.testing_path = self.path / "data" / "fir_conv"

    def test_from_conv_file(self, plot=False):
        """
        Test running FIRConverter.
        """
        # Read in original and fir-converted data
        s_orig = read(str(self.testing_path /
                          "example_orig_LSVEI_20150827061217.mseed"), 'MSEED')
        s_fir = read(str(self.testing_path /
                         "example_FIR_LSVEI_20150827061217.mseed"), 'MSEED')

        # Correct orig data
        obj = FIRConverter.from_conv_file(str(self.testing_path / 'lc2000_fir3_0.json'))
        s_fird, _ = obj.apply(s_orig)
        print(s_orig, s_fir, s_fird)
        s_fird.write('example_FIR_LSVEI_20150827061217.mseed','MSEED')
        if plot is True:
            # Oversample for smooveness
            orig = s_orig[0].copy().resample(500).detrend()
            ref = s_fir[0].copy().resample(500)
            test = s_fird[0].copy().resample(500)
            fig, ax = plt.subplots()
            plt.plot(orig.data, color='k', lw=2, label='original')
            plt.plot(ref.data, color='b', label='reference')
            plt.plot(test.data, color='r', ls='--', label='test')
            plt.legend()
            plt.show()
        self.assertListEqual(s_fir[0].data.tolist(), s_fird[0].data.tolist())

    def test_from_builtin(self):
        """
        Test running FIRConverter.
        """
        # Read in original and fir-converted data
        s_orig = read(str(self.testing_path /
                          "example_orig_LSVEI_20150827061217.mseed"), 'MSEED')
        s_fir = read(str(self.testing_path /
                         "example_FIR_LSVEI_20150827061217.mseed"), 'MSEED')

        s_fird, _ = FIRConverter.from_builtin('lc2000_fir3_0').apply(s_orig)
        self.assertListEqual(s_fir[0].data.tolist(), s_fird[0].data.tolist())

    def test_root_classification(self):
        lin_fir = [0.25, -0.5, 1, -0.5, 0.25]
        fc = FIRConverter.from_zeros(lin_fir, sps=100, timetag=0, decimation_factor=2)
        z = fc.z
        assert len(z) == 4
        assert len(fc.z_max) + len(fc.z_min) + len(fc.z_UC) == 4

    # def test_minimum_phase_filter_normalization():
    #     fir = [1, -0.5, 0.25]
    #     fc = FIRConverter.from_zeros(fir, sps=100, timetag=0)
    #     fir_mp = fc.compute_minimum_phase_filter()
    # 
    #     # Minimum-phase filter should sum to 1 after normalization
    #     assert np.isclose(np.sum(fir_mp), 1.0)

    def test_AR_MA_shapes(self):
        lin_fir = [0.25, -0.5, 1, -0.5, 0.25]
        fc = FIRConverter.from_zeros(lin_fir, sps=100, timetag=0, decimation_factor=2)
        a, b = fc._compute_AR_MA(fc.z)

        assert len(b) == len(a) + 1
        assert np.isfinite(a).all()
        assert np.isfinite(b).all()

def suite():
    return unittest.makeSuite(TestFirConvMethods, 'test')


if __name__ == '__main__':
    unittest.main(defaultTest='suite')
