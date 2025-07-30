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
from copy import deepcopy
from difflib import unified_diff

from obspy.core.inventory.response import FIRResponseStage
from obspy import read_inventory
from obspy.core.stream import read as stream_read
from matplotlib import pyplot as plt
import numpy as np

from tiskitpy import Compliance
from tiskitpy.compliance import EarthModel1D


class TestMethods(unittest.TestCase):
    """
    Test suite
    """
    def setUp(self):
        self.path = Path(inspect.getfile(
            inspect.currentframe())).resolve().parent
        self.test_path = self.path / "data" / "compliance"
        freqs = np.arange(0.005, 0.01, 0.001) # Just a few values for testing
        self.compliance = Compliance(freqs,
                                     np.ones(freqs.shape),
                                     1.e-1*np.ones(freqs.shape),
                                     2300., "output", False)

    def _compare_to_ref_file(self, a, b):
        "a: test_file, b: reference_file"
        with open(b, "r") as f:
            expected_lines = f.readlines()
        with open(a, "r") as f:
            actual_lines = f.readlines()

        print(a)
        print(expected_lines)
        print(b)
        print(actual_lines)
        diff = list(unified_diff(expected_lines, actual_lines))
        assert diff == [], "Unexpected file contents:\n" + "".join(diff)

    def test_str(self):
        self.assertEqual(self.compliance.__str__(),
                         'Compliance object:\n'
                         '  5 frequencies, from 0.005 to 0.009000000000000001 Hz\n'
                         "  water_depth='2300.0'\n"
                         '  noise_channel=output\n'
                         '  gravity_corrected=False')
 
    def test_write(self):
        head = 'test'
        
        self.compliance.write(head)
        fname = head + '_compliance_Pa-1.csv'
        self._compare_to_ref_file(self.path / fname, self.test_path / fname)
        (self.path / fname).unlink()

        self.compliance.write('test', units='m/Pa')
        fname = head + '_compliance_m.Pa-1.csv'
        self._compare_to_ref_file(self.path / fname, self.test_path / fname)
        (self.path / fname).unlink()

        self.compliance.write('test', units='m/s/Pa')
        fname = head + '_compliance_m.s-1.Pa-1.csv'
        self._compare_to_ref_file(self.path / fname, self.test_path / fname)
        (self.path / fname).unlink()

        self.compliance.write('test', units='m/s^2/Pa')
        fname = head + '_compliance_m.s-2.Pa-1.csv'
        self._compare_to_ref_file(self.path / fname, self.test_path / fname)
        (self.path / fname).unlink()

        with self.assertRaises(ValueError):
            self.compliance.write('test', units='haha')

    def test_correct_gravity_terms(self):
        return

    def test_convert_compliance(self):
        c, u = self.compliance._convert_compliance('1/Pa')
        self.assertEqual(list(c), list(self.compliance.values))
        self.assertEqual(list(u), list(self.compliance.uncertainties))
        c_as_mPa = np.array([4596.3589923303625, 3762.2939138026377,
                             3155.9357563870667, 2692.0314906653507,
                             2323.3536020140164])
        om = 2 * np.pi * self.compliance.freqs
        c, u = self.compliance._convert_compliance('m/Pa')
        self.assertEqual(list(c.astype('int')), list((c_as_mPa).astype('int')))
        self.assertEqual(list(u.astype('int')), list((c_as_mPa/10.).astype('int')))
        c, u = self.compliance._convert_compliance('m/s/Pa')
        self.assertEqual(list(c.astype('int')), list((c_as_mPa*om).astype('int')))
        self.assertEqual(list(u.astype('int')), list((c_as_mPa*om/10.).astype('int')))
        c, u = self.compliance._convert_compliance('m/s^2/Pa')
        self.assertEqual(list(c.astype('int')), list((c_as_mPa*om**2).astype('int')))
        self.assertEqual(list(u.astype('int')), list((c_as_mPa*om**2/10.).astype('int')))
        with self.assertRaises(ValueError):
            c, u = self.compliance._convert_compliance('Pa')

    def test_gravd(self):
        gravd = Compliance.gravd
        g = 9.81
        # Shallow water cases
        om = 0.001
        for H in (200., 400.):
            print(gravd([om], H))
            self.assertAlmostEqual(gravd([om], H)[0][0], om/np.sqrt(g*H))
        # Deep water cases
        H = 5000.
        for om in (1., 10., 100.):
            print(f'{om=}')
            self.assertAlmostEqual(gravd([om], H)[0][0], (om**2)/g)
 
    def test_calc_norm_compliance(self):
        """Also tests simple EarthModel1D class"""
        rho=3000
        vp = 6000
        vs = vp/np.sqrt(3)  # Poisson solid
        freqs = np.array([0.001, 0.003, 0.005])
        hs_model = EarthModel1D([[1000, rho, vp, vs],
                                 [1000, rho, vp, vs]])
        theo_norm_compl = - vp**2 / (2 * rho * vs**2 * (vp**2 - vs**2))
        delta = -theo_norm_compl/400  # Require < 0.25% difference
        print(f'{theo_norm_compl=}')
        for H in (10., 100., 1000., 2000., 4000.): 
            # Differences are bigger as water is deeper (ocean waves are faster)
            print(f'{H=}')
            nc = Compliance.calc_norm_compliance(H, freqs, hs_model)
            for x in nc:
                print(f'{np.abs(100*(x-theo_norm_compl)/x):.02f}% difference')
                self.assertAlmostEqual(x, theo_norm_compl, delta=delta)
        # self.assertEqual(x, theo_norm_compl)
                             
    def test_zp_to_ncompl(self):
        gravd = Compliance.gravd
        freqs = np.array([.001])
        zp = np.array([1.])
        H = 2000.
        omega = 2 * np.pi * freqs
        self.assertEqual(Compliance._zp_to_ncompl(freqs, zp, H, 'M'),
                         gravd(omega, H)*zp)
        self.assertAlmostEqual(Compliance._zp_to_ncompl(freqs, zp, H, 'M/S')[0],
                         (gravd(omega, H)*zp/omega)[0])
        self.assertAlmostEqual(Compliance._zp_to_ncompl(freqs, zp, H, 'M/S^2')[0],
                         (gravd(omega, H)*zp/omega**2)[0])
                            

def suite():
    return unittest.makeSuite(TestMethods, 'test')



if __name__ == '__main__':
    unittest.main(defaultTest='suite')
