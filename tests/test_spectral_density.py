#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Functions to test spectral functions
"""
import unittest
import inspect
from pathlib import Path

from obspy import read_inventory
from obspy.core.stream import read as obspy_read

from crawtools.spectral import PSD


class TestMethods(unittest.TestCase):
    """
    Test suite
    """
    def setUp(self):
        self.path = Path(inspect.getfile(
            inspect.currentframe())).resolve().parent
        self.test_path = self.path / "data" / "spectral_density"

    def test_PSD_calc(self):
        stream = obspy_read(str(self.test_path / 'XS.S10D..LHZ.mseed'),
                            'MSEED')
        inv = read_inventory(str(self.test_path / 'stations_PILAB_S10dec.xml'))
        psd = PSD.calc(stream.select(component='Z')[0], inv=inv)
        self.assertEqual(psd.freqs[0], 0.0009765625)
        self.assertEqual(psd.units, 'dB ref 1 (m/s^2)^2/Hz')
        self.assertAlmostEqual(psd.data[0], 9.50923118876e-12)


def suite():
    return unittest.makeSuite(TestMethods, 'test')


if __name__ == '__main__':
    unittest.main(defaultTest='suite')
