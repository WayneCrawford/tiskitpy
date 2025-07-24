#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Functions to test the  functions
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

from tiskitpy import SeafloorSynthetic, Compliance, PSDVals, from_DBs, to_DBs


class TestMethods(unittest.TestCase):
    """
    Test suite
    """
    def setUp(self):
        self.path = Path(inspect.getfile(
            inspect.currentframe())).resolve().parent
        self.test_path = self.path / "data" / "decimate"
        self.seafloor_synthetic = SeafloorSynthetic()  # Uses all defaults

    def test_compliance_noise_IG_Pa_seafloor(self):
        """
        Just check that the lowest frequency value is correct
        """
        H = self.seafloor_synthetic.water_depth
        omega_IG = self.seafloor_synthetic.IG_m_seasurface.freqs*2*np.pi
        k = Compliance.gravd(omega_IG, H)
        seawater_density = 1030  #  1020-1029 at the surface, up to 1050 at deep seafloor
        g = 9.81  # 9.78 at equator, 9.83 at poles
        Pa_per_m = seawater_density*g
        self.assertAlmostEqual(self.seafloor_synthetic.IG_Pa_seafloor.values[0],
                               self.seafloor_synthetic.IG_m_seasurface.values[0]
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
