#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Functions to test the FIRConverter class
"""
import unittest
import inspect
import tempfile
from pathlib import Path

from obspy import read, Stream, Trace, UTCDateTime
from obspy.core.inventory.response import FIRResponseStage
from scipy.signal import freqz
from matplotlib import pyplot as plt
import numpy as np

from tiskitpy import FIRConverter


class TestFirConvMethods(unittest.TestCase):
    """
    Test suite for FIRConverter operations.
    """
    def setUp(self):
        self.path = Path(inspect.getfile(
            inspect.currentframe())).resolve().parent
        print(self.path)
        self.testing_path = self.path / "data" / "fir_conv"
        self.stream_orig = read(
            str(self.testing_path / "example_orig_LSVEI_20150827061217.mseed"),
            'MSEED')
        self.stream_fir = read(
            str(self.testing_path / "example_FIR_LSVEI_20150827061217.mseed"),
            'MSEED')
        self.conv_file = str(self.testing_path / 'lc2000_fir3_0.json')
        self.builtin = 'lc2000_fir3_0'

    def test_from_conv_file(self, plot=False):
        """
        Test running FIRConverter.
        """
        # Read in original and fir-converted data
        converter = FIRConverter.from_conv_file(self.conv_file)
        stream_fir, _ = converter.apply(self.stream_orig)
        print(self.stream_orig, self.stream_fir, stream_fir)
        self.assertListEqual(self.stream_fir[0].data.tolist(),
                             stream_fir[0].data.tolist())
        stream_fir.write('example_FIR_LSVEI_20150827061217.mseed','MSEED')
        if plot is True:
            self.plot_oversampled([
                [self.stream_orig[0], {'color': 'k', 'lw': 2, 'label': 'original'}],
                [self.stream_fir[0], {'color': 'k', 'label': 'reference'}],
                [stream_fir[0], {'color': 'r', 'ls': '--', 'label': 'test'}]
                             ])

    def plot_oversampled(self, trace_lists):
            fig, ax = plt.subplots()
            for tl in trace_lists:
                # Oversample for smooveness
                trace = tl[0].copy()
                trace = trace.resample(10*trace.stats.sampling_rate)
                kwargs = tl[1]
                plt.plot(trace.detrend(), **kwargs)
            plt.legend()
            plt.show()

    def test_apply_buffer(self, plot=True):
        # Test applying the code to a local EQ with visible FIR effect
        # Compare working on the entire stream versus the entire stream
        # chopped in two
        
        # Cut the sample data in two
        assert len(self.stream_orig) == 1 # Make sure only 1, contiguous trace
        t = self.stream_orig[0]
        sr = t.stats.sampling_rate
        midtime = t.stats.starttime + (t.stats.endtime-t.stats.starttime)/2
        stream_a = self.stream_orig.copy().slice(None, midtime,
                                                 nearest_sample=False)
        stream_b = self.stream_orig.copy().slice(midtime+1/(100*sr), None,
                                                 nearest_sample=False)
        
        # Convert the two halves, stitch them together, and confirm all is good
        converter = FIRConverter.from_conv_file(self.conv_file)
        stream_a_fir, buffer = converter.apply(stream_a, next_stream=stream_b,
                                               detrend=False)
        stream_b_fir, buffer = converter.apply(stream_b, buffer, detrend=False)
        stream_fir, _ = converter.apply(self.stream_orig, detrend=False)
        stream_ab_fir = stream_a_fir + stream_b_fir
        stream_ab_fir.merge()
        print(stream_ab_fir)
        if plot is True:
            self.plot_oversampled([
                [self.stream_orig[0], {'color': 'k', 'lw': 2, 'label': 'original'}],
                [self.stream_fir[0], {'color': 'g', 'label': 'reference'}],
                [stream_ab_fir[0], {'color': 'r', 'ls': '--', 'label': 'test_cut'}],
                [stream_fir[0], {'color': 'b', 'ls': '--', 'label': 'test_uncut'}]
                             ])
        self.assertListEqual(stream_fir[0].data.tolist(),
                             stream_ab_fir[0].data.tolist())
        # Output has the same samples as the input
        self.assertEqual(stream_a_fir[0].stats.starttime,
                         stream_a[0].stats.starttime)
        self.assertEqual(stream_a_fir[0].stats.npts, stream_a[0].stats.npts)

        # Default (detrend=None) only detrends traces without neighbours
        stream_a_auto, _ = converter.apply(stream_a, next_stream=stream_b)
        self.assertListEqual(stream_a_fir[0].data.tolist(),
                             stream_a_auto[0].data.tolist())
        self.assertListEqual(
            converter.apply(stream_a, detrend=True)[0][0].data.tolist(),
            converter.apply(stream_a)[0][0].data.tolist())

        # Second half shifted off the sample grid (e.g. clock correction):
        # by default its samples are treated as consecutive
        stream_b_shift = stream_b.copy()
        stream_b_shift[0].stats.starttime += 0.06 / sr
        stream_a_fir_shift, buffer = converter.apply(
            stream_a, next_stream=stream_b_shift, detrend=False)
        stream_b_fir_shift, _ = converter.apply(stream_b_shift, buffer,
                                                detrend=False)
        self.assertListEqual(
            stream_fir[0].data.tolist(),
            stream_a_fir_shift[0].data.tolist()
            + stream_b_fir_shift[0].data.tolist())
        # unless max_misalignment is smaller than the shift
        stream_b_fir_shift, _ = converter.apply(
            stream_b_shift, buffer, detrend=False, max_misalignment=0.01)
        self.assertNotEqual(stream_b_fir[0].data.tolist(),
                            stream_b_fir_shift[0].data.tolist())

        # Both halves in one stream
        stream_ab_fir, _ = converter.apply(stream_a + stream_b, detrend=False)
        self.assertListEqual(stream_fir[0].data.tolist(),
                             stream_ab_fir.merge()[0].data.tolist())

        # Make a gap in the data, confirm that it is correctly handled
        stream_a[0].data = stream_a[0].data[0:-5]
        stream_b_alone, _ = converter.apply(stream_b, detrend=False)
        for gap_handling in ('interpolate', 'edge'):
            stream_a_fir, buffer = converter.apply(
                stream_a, next_stream=stream_b, detrend=False,
                gap_handling=gap_handling)
            stream_b_fir, _ = converter.apply(
                stream_b, buffer, detrend=False, gap_handling=gap_handling)
            self.assertEqual(stream_a_fir[0].stats.endtime,
                             stream_a[0].stats.endtime)
            self.assertEqual(stream_b_fir[0].stats.starttime,
                             stream_b[0].stats.starttime)
            if gap_handling == 'edge':
                self.assertListEqual(stream_b_alone[0].data.tolist(),
                                     stream_b_fir[0].data.tolist())
            else:
                self.assertNotEqual(stream_b_alone[0].data.tolist(),
                                    stream_b_fir[0].data.tolist())
            if plot is True:
                self.plot_oversampled([
                    [self.stream_orig[0], {'color': 'k', 'lw': 2, 'label': 'original'}],
                    [self.stream_fir[0], {'color': 'k', 'label': 'reference'}],
                    [stream_a_fir[0], {'color': 'r', 'ls': '--', 'label': 'test_a'}],
                    [stream_b_fir[0], {'color': 'b', 'ls': '--', 'label': 'test_b'}]
                                 ])

        # Gap longer than max_interp_gap is not interpolated
        stream_b_fir, _ = converter.apply(stream_b, stream_a, detrend=False,
                                          max_interp_gap=4)
        self.assertListEqual(stream_b_alone[0].data.tolist(),
                             stream_b_fir[0].data.tolist())

    def test_apply_SDS(self):
        # Make an SDS archive of day files crossing a year boundary
        rng = np.random.default_rng(42)
        sr = 1.
        t0 = UTCDateTime('2015-12-30')
        tr = Trace(rng.integers(-10000, 10000, int(4*86400*sr), dtype='int32'),
                   header={'network': 'XX', 'station': 'STA',
                           'location': '00', 'channel': 'BHZ',
                           'sampling_rate': sr, 'starttime': t0})
        converter = FIRConverter.from_conv_file(self.conv_file)
        with tempfile.TemporaryDirectory() as tmpdir:
            sds_in = Path(tmpdir) / 'SDS_in'
            sds_out = Path(tmpdir) / 'SDS_out'
            for day in range(4):
                t = t0 + day*86400
                day_tr = tr.slice(t, t + 86400 - 1/sr)
                fpath = (sds_in / str(t.year) / 'XX' / 'STA' / 'BHZ.D' /
                         f'XX.STA.00.BHZ.D.{t.year}.{t.julday:03d}')
                fpath.parent.mkdir(parents=True, exist_ok=True)
                day_tr.write(str(fpath), format='MSEED')
            # Skip the last day, which should still be used for the 3rd day
            converter.apply_SDS(sds_in, sds_out,
                                enddate=t0 + 3*86400, loc_code='01')
            out_files = sorted(x.name for x in sds_out.glob('*/*/*/*/*'))
            self.assertListEqual(out_files, ['XX.STA.01.BHZ.D.2015.364',
                                             'XX.STA.01.BHZ.D.2015.365',
                                             'XX.STA.01.BHZ.D.2016.001'])
            st_out = read(str(sds_out / '*/*/*/*/*')).merge()
            self.assertEqual(len(st_out), 1)
            self.assertEqual(st_out[0].id, 'XX.STA.01.BHZ')
            # Same as converting the whole trace in one piece
            st_whole, _ = converter.apply(Stream([tr]), detrend=False)
            st_whole.trim(t0, t0 + 3*86400 - 1/sr)
            self.assertListEqual(st_whole[0].data.tolist(),
                                 st_out[0].data.tolist())
            # loc_code already in output directory
            with self.assertRaises(ValueError):
                converter.apply_SDS(sds_in, sds_out, loc_code='01')

    def test_response_stage(self):
        converter = FIRConverter.from_conv_file(self.conv_file)
        values = converter.output_stage_values()
        self.assertEqual(values['symmetry'], 'NONE')
        self.assertEqual(values['delay'], 0.)
        self.assertEqual(values['correction'], 0.)
        self.assertAlmostEqual(sum(values['coefficients']), 1.)
        # Linear phase input stage, given as half of an ODD symmetric FIR
        # with integer coefficients
        fs_fir = 125.
        firzeros = converter.firzeros
        delay = converter.timetag / fs_fir
        stage = FIRResponseStage(
            3, 1., 1., 'COUNTS', 'COUNTS', symmetry='ODD', name='FIR3',
            coefficients=firzeros[:len(firzeros)//2 + 1].tolist(),
            decimation_input_sample_rate=fs_fir, decimation_factor=2,
            decimation_offset=0, decimation_delay=delay,
            decimation_correction=delay)
        new_stage = converter.convert_response_stage(stage)
        self.assertEqual(new_stage.symmetry, 'NONE')
        self.assertEqual(new_stage.decimation_delay, 0.)
        self.assertEqual(new_stage.decimation_correction, 0.)
        self.assertEqual(new_stage.name, 'FIR3 (minimum phase)')
        self.assertEqual(stage.symmetry, 'ODD')  # input stage unchanged
        new_coeffs = np.array([float(x) for x in new_stage.coefficients])
        half = firzeros[:len(firzeros)//2 + 1]
        self.assertAlmostEqual(new_coeffs.sum() / (2*half.sum() - half[-1]),
                               1.)
        # Same amplitude response, minimum phase
        _, h_in = freqz(firzeros, worN=512)
        _, h_out = freqz(new_coeffs, worN=512)
        np.testing.assert_allclose(np.abs(h_out), np.abs(h_in),
                                   atol=1e-6 * np.abs(h_in).max())
        self.assertTrue(np.all(np.abs(np.roots(new_coeffs)) <= 1.001))
        # Wrong FIR
        stage.coefficients = stage.coefficients[:-1]
        with self.assertRaises(ValueError):
            converter.convert_response_stage(stage)

    def test_pad_lengths(self):
        converter = FIRConverter.from_conv_file(self.conv_file)
        # 25 samples of timetag shift + 11 of interpolation filter
        self.assertEqual(converter.pad_lengths()[0], 36)
        n_after = converter.pad_lengths()[1]
        self.assertGreater(converter.pad_lengths(rel_tol=1e-15)[1], n_after)
        self.assertLess(converter.pad_lengths(rel_tol=1e-3)[1], n_after)


    def test_from_builtin(self):
        """
        Test running FIRConverter.
        """
        # Read in original and fir-converted data
        converter = FIRConverter.from_builtin(self.builtin)
        stream_fir, _ = converter.apply(self.stream_orig)
        self.assertListEqual(self.stream_fir[0].data.tolist(),
                             stream_fir[0].data.tolist())

    def test_root_classification(self):
        lin_fir = [0.25, -0.5, 1, -0.5, 0.25]
        fc = FIRConverter.from_zeros(lin_fir, sps=100, timetag=0, decimation_factor=2)
        z = fc.z
        assert len(z) == 4
        assert len(fc.z_max) + len(fc.z_min) + len(fc.z_UC) == 4


    # def test_minimum_phase_filter_normalization():
    #     fir = [1, -0.5, 0.25]
    #     fc = FIRConverter.from_zeros(fir, sps=100, timetag=0)
    #     fir_mp = fc.compute_minimum_phase_filter()
    # 
    #     # Minimum-phase filter should sum to 1 after normalization
    #     assert np.isclose(np.sum(fir_mp), 1.0)

    def test_AR_MA_shapes(self):
        lin_fir = [0.25, -0.5, 1, -0.5, 0.25]
        fc = FIRConverter.from_zeros(lin_fir, sps=100, timetag=0, decimation_factor=2)
        a, b = fc._compute_AR_MA(fc.z)

        assert len(b) == len(a) + 1
        assert np.isfinite(a).all()
        assert np.isfinite(b).all()

def suite():
    return unittest.makeSuite(TestFirConvMethods, 'test')


if __name__ == '__main__':
    unittest.main(defaultTest='suite')
