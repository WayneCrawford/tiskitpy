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

from tiskitpy import CleanRotator, CleanSequence as CS


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
        rotator = CleanRotator(
            stream, verbose=False,
            remove_eq=str(self.test_path /
                          "20161205T-20161207T_MM5.85_eqcat.qml"),
            save_eq_file=False)
        self.assertAlmostEqual(rotator.angle, 0.18, delta=0.01)
        self.assertAlmostEqual(rotator.azimuth, 241.67, delta=0.01)
        
        rot_stream = rotator.apply(stream)
        rotZ = rot_stream.select(channel='*Z')[0]
        self.assertEqual(rotZ.get_id(), 'XS.S10D..LHZ')
        self.assertEqual(CS.seedid_tag(rotZ).get_id(), 'XS.S10D.-ROT.LHZ')
        self.assertEqual(CS.seedid_tag(rotZ, 'min_code').get_id(), 'XS.S10D.-ROT.LHZ')


def suite():
    return unittest.makeSuite(TestMethods, 'test')


if __name__ == '__main__':
    unittest.main(defaultTest='suite')
