#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Functions to test the lcheapo functions
"""
# from os import system
import unittest
import inspect
from pathlib import Path

from obspy.core.stream import read as stream_read

from tiskit import CleanRotator


class TestMethods(unittest.TestCase):
    """
    Test suite
    """
    def setUp(self):
        self.path = Path(inspect.getfile(
            inspect.currentframe())).resolve().parent
        self.test_path = self.path / "data" / "rotate_clean"

    def test_clean_rotator(self):
        """
        Test CleanRotator class.
        """
        # this file only has 3 minutes and an EQ, is it enough?
        stream = stream_read(str(self.test_path / 'XS.S10D.LH.mseed'))
        # rot_stream, bestAngle, bestAzimuth = rotate_clean(stream,
        #                                                   verbose=False)
        rotator = CleanRotator(stream, verbose=False)
        # rot_stream = rotator.apply(stream)
        self.assertAlmostEqual(rotator.angle, -0.18, delta=0.01)
        self.assertAlmostEqual(rotator.azimuth, 61.67, delta=0.01)


def suite():
    return unittest.makeSuite(TestMethods, 'test')


if __name__ == '__main__':
    unittest.main(defaultTest='suite')
