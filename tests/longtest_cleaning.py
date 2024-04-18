#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Functions to test spectral functions
"""
import unittest
from pathlib import Path
import inspect
import pickle

from obspy.core.stream import read
from obspy.core.inventory import read_inventory
import numpy as np
from matplotlib import pyplot as plt

from tiskitpy import SpectralDensity, DataCleaner, CleanRotator
from make_test_stream import make_test_stream


class TestMethods(unittest.TestCase):
    """
    Test suite
    """

    def setUp(self):
        self.path = (Path(inspect.getfile(inspect.currentframe()))
                     .resolve().parent)
        self.test_path = self.path / "data" / "cleaning"

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
        rotator = CleanRotator(
            stream, filt_band=rotate_band, verbose=False, uselogvar=False,
            remove_eq=str(self.test_path / "test_eqcat.qml"))
        stream_rot = rotator.apply(stream)
        dc = DataCleaner(stream, dc_tforder, **dc_kwargs)
        dc_rot = DataCleaner(stream_rot, dc_tforder, **dc_kwargs)
        stream_dced = dc.clean_stream(stream, in_time_domain=True)
        stream_rot_dced = dc_rot.clean_stream(stream_rot, in_time_domain=True)
        fig, ax = plt.subplots(1)
        psds=[]
        for stream, label, minval in zip(
                (stream, stream_rot, stream_dced, stream_rot_dced),
                ('original', 'rotated', 'stream_corr', 'rot + stream_corr'),
                (-164.2, -166.5, -166.1, -170.7)):
            st = stream.select(component='Z')
            st[0].stats.channel='LHZ'
            st[0].stats.location=''
            sd=SpectralDensity.from_stream(st, inv=inv, window_s=2048)
            psds.append(10*np.log10(sd.autospect('XS.S11D..LHZ')))
            self.assertAlmostEqual(np.min(psds[-1]), minval, delta=0.1)
            ax.semilogx(sd.freqs, psds[-1], label=label)
        # test SpectralDensity calculation with datacorrection provided
        sd=SpectralDensity.from_stream(stream_rot, window_s=2048, inv=inv,
                                       data_cleaner=dc_rot)
        psds.append(10*np.log10(sd.autospect('XS.S11D..LHZ-1-2-H')))
        self.assertAlmostEqual(np.min(psds[-1]), -175.5, delta=0.1)
        ax.semilogx(sd.freqs, psds[-1], label='rot + psd_corr')
        ax.legend()
        # plt.savefig('longtest_cleaning.png')
        # with open('psds.p', 'wb') as o:
        #     pickle.dump(psds, o)
        with open(str(self.test_path / 'psds.p'), 'rb') as i:
            ref_psds = pickle.load(i)
        for psd, ref_psd in zip(psds, ref_psds):
            self.assertTrue(np.all(psd == ref_psd))
        
        

def suite():
    return unittest.makeSuite(TestMethods, "test")


if __name__ == "__main__":
    unittest.main(defaultTest="suite")
