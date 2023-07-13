#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Functions to test spectral functions
"""
import unittest
# from pathlib import Path
# import inspect

from obspy.core.stream import read

from tiskitpy import CleanerString as CS


class TestMethods(unittest.TestCase):
    """
    Test suite
    """

    def setUp(self):
        # self.path = (Path(inspect.getfile(inspect.currentframe()))
        #              .resolve().parent)
        # self.test_path = self.path / "data" / "data_cleaner"
        self.window_s = 100.0

    def test_make(self):
        """Test make()"""
        self.assertEqual(CS.make('BHZ', '-X-Y'), '-Z')
        self.assertEqual(CS.make('BHX', '-X-Y'), '-HX')
        with self.assertRaises(ValueError):
            CS.make('X', '-X-Y')

    def test_insert(self):
        """Test insert()"""
        self.assertEqual(CS.insert('-X-Y', 'BHZ'), '-X-Y.BHZ')
        self.assertEqual(CS.insert('-X-Y', 'XX.STA.00.BHZ'), 'XX.STA.00-X-Y.BHZ')
        
    def test_strip(self):
        """Test strip()"""
        self.assertEqual(CS.strip('-X-Y.BHZ'), 'BHZ')
        self.assertEqual(CS.strip('00-X-Y.BHZ'), '00.BHZ')
        self.assertEqual(CS.strip('.BHZ'), '.BHZ')
        self.assertEqual(CS.strip('XX.STA.00-X-Y.BHZ'), 'XX.STA.00.BHZ')
        self.assertEqual(CS.strip('XX.STA.-X-Y.BHZ'), 'XX.STA..BHZ')
        self.assertEqual(CS.strip('XX.STA.00.BHZ'), 'XX.STA.00.BHZ')
        

def suite():
    return unittest.makeSuite(TestMethods, "test")


if __name__ == "__main__":
    unittest.main(defaultTest="suite")
