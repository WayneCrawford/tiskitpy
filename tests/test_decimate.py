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

from obspy.core.inventory.response import FIRResponseStage
from obspy import read_inventory
from obspy.core.stream import read as stream_read
from matplotlib import pyplot as plt
import numpy as np

from tiskitpy import Decimator
from tiskitpy.decimate import FIRFilter


class TestMethods(unittest.TestCase):
    """
    Test suite
    """
    def setUp(self):
        self.path = Path(inspect.getfile(
            inspect.currentframe())).resolve().parent
        self.test_path = self.path / "data" / "decimate"

    def assertTextFilesEqual(self, first, second, msg=None):
        with open(first) as f:
            str_a = f.read()
        with open(second) as f:
            str_b = f.read()

        if str_a != str_b:
            first_lines = str_a.splitlines(True)
            second_lines = str_b.splitlines(True)
            delta = difflib.unified_diff(
                first_lines, second_lines,
                fromfile=first, tofile=second)
            message = ''.join(delta)

            if msg:
                message += " : " + msg

            self.fail("Multi-line strings are unequal:\n" + message)

    def assertBinFilesEqual(self, first, second, msg=None):
        """ Compare two binary files """
        self.assertTrue(filecmp.cmp(first, second))

    def test_FIR_Filter(self):
        """Test FIRFilter outputs"""
        coeffs = [0, -0.5, 2, -0.5, 0]
        decim = 2
        offset = 2
        name = 'fake'
        sampling_rate = 100
        ff = FIRFilter(coeffs, offset=offset, decim=decim, name=name)
        ff_obspy = ff.to_obspy(sampling_rate)
        ff_test = FIRResponseStage(
            stage_sequence_number=0,
            stage_gain=1,
            stage_gain_frequency=0,
            input_units='count',
            input_units_description='digital counts',
            output_units='count',
            output_units_description='digital counts',
            symmetry='NONE',
            name=name,
            coefficients=coeffs,
            decimation_input_sample_rate=sampling_rate,
            decimation_factor=decim,
            decimation_offset=offset,
            decimation_delay=offset/sampling_rate,
            decimation_correction=offset/sampling_rate)
        # print(ff_obspy)
        # print(ff_test)
        self.assertEqual(ff_test, ff_obspy)

    def test_plot_SAC_FIR_Filters(self):
        """ Test"""
        plot = False
        for decim in range(2, 8):
            ff = FIRFilter.from_SAC(decim)
            if plot:
                ff.plot()
                plt.savefig(f'SAC_decim{decim}.png')

    def test_decimate_data(self):
        """Run decimation on data with a big earthquake"""
        st = stream_read(str(self.test_path / 'XS.S10D.20161212T2053.mseed'))
        decim = Decimator([5])
        st2 = decim.decimate(st)
        self.assertTrue(len(st[0].data == 5 * len(st2[0].data)))

        # Test dtype handling
        self.assertEqual(st2[0].data.dtype, st[0].data.dtype )
        st3 = decim.decimate(st, keep_dtype=False)
        self.assertNotEqual(st3[0].data.dtype, st[0].data.dtype )

    def test_decimate_impulse(self):
        """Run decimation on data with a single impulse"""
        plot = False
        for d in range(2, 8):
            st = stream_read(str(self.test_path
                                 / 'XS.S10D.20161212T2053.mseed'))
            tr = st.select(channel='BHZ')[0]
            tr.data = np.zeros(tr.data.shape)
            tr.data[int(len(tr.data)/2)] = 1
            decim = Decimator([d])
            tr2 = decim.decimate(tr)
            if plot:
                plot_compare(tr, tr2, d,
                             '2016-12-12T20:54:27', '2016-12-12T20:54:33',
                             savefig=f'test_SACdecim{d}_impulse.png')

    def test_SAC_Filters(self):
        """ Test lcdump outputs"""
        for decim in range(2, 8):
            ff = FIRFilter.from_SAC(decim, normalize=False)
            fname = f'decim{decim}.pickle'
            # The following is to save a file to put in self.test_path
            with open(fname, 'wb') as f:
                pickle.dump(ff, f)
            with open(str(self.test_path / fname), 'rb') as f:
                fref = pickle.load(f)
            self.assertEqual(ff, fref)
            Path(fname).unlink()

    def test_bad_SAC(self):
        """
        Should give error if non-existant decimation demanded
        """
        args = [8]
        self.assertRaises(ValueError, FIRFilter.from_SAC, *args)

    def test_modify_stationxml(self):
        """Read a stationxml file and create new decimated data channels"""
        inv = read_inventory(str(self.test_path / 'stations_PILAB.xml'))
        # Restrict to Scripps stations to shorten printouts
        inv_SIO = inv.select(station='S1*D')
        decim = Decimator([2, 5, 5])
        st = stream_read(str(self.test_path / 'XS.S10D.20161212T2053.mseed'))
        # Test using stream channels
        # SIO's FIR8 causes problems here, because it's sum[coeffs] = 0.9767
        invnew = decim.update_inventory(inv_SIO, st, quiet=True)
        self.assertNotEqual(inv_SIO, invnew)
        invnew.write('stations_PILAB_S10decim.xml', format='STATIONXML')
        # self.assertEqual(
        #     inv_SIO, read_inventory(str(self.test_path /
        #                                 'stations_PILAB_S10decim.xml')))
        # Test using explicit channels
        invnew = decim.update_inventory_from_nslc(inv_SIO, 'S10D', 'BHZ',
                                                  quiet=True)
        self.assertNotEqual(inv, invnew)
        invnew.write('stations_PILAB_S10BHZdecim.xml', format='STATIONXML')
        # self.assertEqual(
        #     inv_SIO, read_inventory(str(self.test_path /
        #                                 'stations_PILAB_S10decim.xml')))
        Path('stations_PILAB_S10decim.xml').unlink()
        Path('stations_PILAB_S10BHZdecim.xml').unlink()


def suite():
    return unittest.makeSuite(TestMethods, 'test')


def plot_compare(tr1, tr2, decim, startdate, enddate, show=False,
                 savefig='plot_compare.png'):
    """
    Plot a comparison of pre- and post-decimation_correction

    Args:
        tr1 (:class:`obspy.core.stream.Trace`): first trace
        tr2 (:class:`obspy.core.stream.Trace`): second trace
        startdate (str): in datetime isoformat
        enddate (str): in datetime isoformat
        show (bool): plot on screen
        savefig (str): filename to save plot to
    """
    # print(f'saving to {savefig}')
    plt.plot_date(tr1.times(type='matplotlib'), tr1.data, 'r-',
                  label='original')
    plt.plot(tr2.times(type='matplotlib'), tr2.data, 'b-',
             label=f'decim{decim}')
    plt.xlim(datetime.datetime.fromisoformat(startdate),
             datetime.datetime.fromisoformat(enddate))
    plt.legend()
    if show:
        plt.show()
    if savefig is not None:
        plt.savefig(savefig)
    plt.close()


if __name__ == '__main__':
    unittest.main(defaultTest='suite')
