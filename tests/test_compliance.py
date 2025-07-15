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

from obspy.core.inventory.response import FIRResponseStage
from obspy import read_inventory
from obspy.core.stream import read as stream_read
from matplotlib import pyplot as plt
import numpy as np

from tiskitpy.compliance import (ComplianceNoise, PSDVals, EarthModel1D, 
                                 gravd, calc_norm_compliance,
                                 zp_to_norm_compliance, from_DBs, to_DBs)


class TestMethods(unittest.TestCase):
    """
    Test suite
    """
    def setUp(self):
        self.path = Path(inspect.getfile(
            inspect.currentframe())).resolve().parent
        self.test_path = self.path / "data" / "decimate"
        self.compliance_noise = ComplianceNoise()  # Uses all defaults

    def test_gravd(self):
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
            nc = calc_norm_compliance(H, freqs, hs_model)
            for x in nc:
                print(f'{np.abs(100*(x-theo_norm_compl)/x):.02f}% difference')
                self.assertAlmostEqual(x, theo_norm_compl, delta=delta)
        # self.assertEqual(x, theo_norm_compl)
                             
    def test_zp_to_norm_compliance(self):
        freqs = np.array([.001])
        zp = np.array([1.])
        H = 2000.
        omega = 2 * np.pi * freqs
        self.assertEqual(zp_to_norm_compliance(freqs, zp, H, 'M'),
                         gravd(omega, H)*zp)
        self.assertAlmostEqual(zp_to_norm_compliance(freqs, zp, H, 'M/S')[0],
                         (gravd(omega, H)*zp/omega)[0])
        self.assertAlmostEqual(zp_to_norm_compliance(freqs, zp, H, 'M/S^2')[0],
                         (gravd(omega, H)*zp/omega**2)[0])
                            
    def test_compliance_noise_IG_Pa_seafloor(self):
        """
        Just check that the lowest frequency value is correct
        """
        H = self.compliance_noise.water_depth
        omega_IG = self.compliance_noise.IG_m_seasurface.freqs*2*np.pi
        k = gravd(omega_IG, H)
        seawater_density = 1030  #  1020-1029 at the surface, up to 1050 at deep seafloor
        g = 9.81  # 9.78 at equator, 9.83 at poles
        Pa_per_m = seawater_density*g
        self.assertAlmostEqual(self.compliance_noise.IG_Pa_seafloor.values[0],
                               self.compliance_noise.IG_m_seasurface.values[0]
                               + 20*np.log10(Pa_per_m/np.cosh(k*H))[0])

    def test_to_DBs(self):
        self.assertEqual(to_DBs(10.),  20.)
        self.assertEqual(to_DBs(100.), 40.)
        self.assertEqual(to_DBs(0.1), -20.)

    def test_from_DBs(self):
        self.assertEqual(from_DBs( 20.), 10.)
        self.assertEqual(from_DBs( 40.), 100.)
        self.assertEqual(from_DBs(-20.), 0.1)
        for x in (0.03, 0.4, 4.555, np.pi, 246, 1.23e5):
            self.assertAlmostEqual(from_DBs(to_DBs(x)), x, delta=x/1e8)

def suite():
    return unittest.makeSuite(TestMethods, 'test')



if __name__ == '__main__':
    unittest.main(defaultTest='suite')
