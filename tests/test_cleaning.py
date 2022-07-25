#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Functions to test spectral functions
"""
import unittest
from pathlib import Path
import inspect

from obspy.core.stream import read
from obspy.core.inventory import read_inventory
import numpy as np

from tiskit import SpectralDensity, DataCleaner, CleanRotator
from make_test_stream import make_test_stream


class TestMethods(unittest.TestCase):
    """
    Test suite
    """

    def setUp(self):
        self.path = (Path(inspect.getfile(inspect.currentframe()))
                     .resolve().parent)
        self.test_path = self.path / "data" / "data_cleaner"
        self.window_s = 100.0

        self.stream, self.sineparms = make_test_stream()
        self.dc = DataCleaner(self.stream, remove_list=["XX.STA.00.BX1"],
                              window_s=self.window_s)

    def test_clean_stream_values(self):
        """Test clean_stream function on real data, compare to reference values"""
        sp_kwargs = {'window_s': 2048, 'windowtype': 'prol1pi'}
        dc_kwargs = {'max_freq': 0.1, 'show_progress': False, **sp_kwargs}
        dc_tforder = ('*1', '*2', '*H')
        rotate_band = (0.003, 0.02)  # freq band in which to calculate rotate angle
        stream = read(str(self.test_path / 'XS.S11D.LH.2016.12.11.mseed'),
                      'MSEED')
        inv = read_inventory(str(self.test_path
                                 / 'stations_PILAB_S_decimated.xml'),
                             'STATIONXML')
        rotator = CleanRotator(stream, filt_band=rotate_band,
                           verbose=False, uselogvar=False)
        stream_rot = rotator.apply(stream)
        dc = DataCleaner(stream, dc_tforder, **dc_kwargs)
        dc_rot = DataCleaner(stream_rot, dc_tforder, **dc_kwargs)
        stream_dced = dc.clean_stream(stream, in_time_domain=True)
        stream_rot_dced = dc_rot.clean_stream(stream_rot, in_time_domain=True)
        for stream, minval in zip((stream, stream_rot, stream_dced, stream_rot_dced),
                                  (-164.2, -166.5, -166.1, -170.7)):
            st = stream.select(component='Z')
            st[0].stats.channel='LHZ'
            st[0].stats.location=''
            sd=SpectralDensity.from_stream(st, inv=inv, window_s=2048)
            min_psd = 10*np.log10(np.min(sd.autospect('XS.S11D..LHZ')))
            self.assertAlmostEqual(min_psd, minval, delta=0.1)
        

def suite():
    return unittest.makeSuite(TestMethods, "test")


if __name__ == "__main__":
    unittest.main(defaultTest="suite")
