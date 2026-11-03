#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Test the SDS console scripts (decimate_SDS, FIR_convert_SDS)
"""
import unittest
import tempfile
from argparse import ArgumentParser
from pathlib import Path

import numpy as np
from obspy import read, read_inventory, Stream, Trace, UTCDateTime
from obspy.core.inventory import (Inventory, Network, Station, Channel, Site,
                                  Response, FIRResponseStage)

from tiskitpy import FIRConverter, Decimator
from tiskitpy.console_scripts import decimate_SDS, fir_convert


class TestSDSScripts(unittest.TestCase):
    """
    Process a synthetic SDS archive whose day files cross a year boundary
    and are offset from each other by 0.06 samples (like clock-corrected
    data), and compare with processing the whole record in one piece
    """
    sr = 1.
    ndays = 3

    def setUp(self):
        self.tmpdir = tempfile.TemporaryDirectory()
        self.root = Path(self.tmpdir.name)
        self.converter = FIRConverter.from_builtin('lc2000_fir3_0')
        rng = np.random.default_rng(3)
        n = int(86400 * self.sr)
        x = np.cumsum(rng.normal(0, 200, self.ndays * n))
        x -= np.convolve(x, np.ones(201) / 201, mode='same')
        data = x.astype('int32')
        t0 = UTCDateTime('2015-12-31T00:00:00.3')
        self.traces = []
        for d in range(self.ndays):
            tr = Trace(data[d*n:(d+1)*n].copy(),
                       header=dict(network='XX', station='STA',
                                   location='00', channel='LHZ',
                                   sampling_rate=self.sr,
                                   starttime=t0 + d*86400 + d*0.06/self.sr))
            t = tr.stats.starttime
            fpath = (self.root / 'SDS' / str(t.year) / 'XX' / 'STA' / 'LHZ.D'
                     / f'XX.STA.00.LHZ.D.{t.year}.{t.julday:03d}')
            fpath.parent.mkdir(parents=True, exist_ok=True)
            tr.write(str(fpath), format='MSEED', encoding='STEIM2')
            self.traces.append(tr)
        self.whole = Trace(
            np.concatenate([t.data for t in self.traces]),
            header={k: self.traces[0].stats[k]
                    for k in ('network', 'station', 'location', 'channel',
                              'starttime', 'sampling_rate')})
        # Inventory: seismometer + LC2000 FIR
        resp = Response.from_paz(zeros=[0j, 0j],
                                 poles=[-0.037+0.037j, -0.037-0.037j],
                                 stage_gain=1500., input_units='M/S',
                                 output_units='V')
        fir = self.converter.firzeros / self.converter.firzeros.sum()
        delay = self.converter.timetag / (2 * self.sr)
        resp.response_stages.append(FIRResponseStage(
            2, 1e6, 0.1, 'V', 'COUNTS', symmetry='ODD', name='FIR3',
            coefficients=fir[:len(fir)//2 + 1].tolist(),
            decimation_input_sample_rate=2*self.sr, decimation_factor=2,
            decimation_offset=0, decimation_delay=delay,
            decimation_correction=delay))
        resp.recalculate_overall_sensitivity(0.1)
        cha = Channel('LHZ', '00', 0., 0., 0., 0., sample_rate=self.sr,
                      response=resp)
        inv = Inventory([Network('XX', [Station('STA', 0., 0., 0.,
                                                channels=[cha],
                                                site=Site('test'))])],
                        source='test')
        self.inv_file = str(self.root / 'stations.xml')
        inv.write(self.inv_file, format='STATIONXML')

    def tearDown(self):
        self.tmpdir.cleanup()

    def _files(self, pattern):
        return [read(str(p))[0] for p in sorted(
            self.root.glob(pattern), key=lambda p: p.name.split('.')[-2:])]

    @staticmethod
    def _args(module, arglist):
        parser = ArgumentParser()
        module.add_arguments(parser)
        return parser.parse_args(arglist)

    def test_FIR_convert_SDS(self):
        fir_convert.main(self._args(fir_convert, [
            str(self.root / 'SDS'), self.inv_file,
            '--builtin', 'lc2000_fir3_0', '-q']))
        out = self._files('SDS/*/XX/STA/LHZ.D/XX.STA.01.*')
        self.assertEqual(len(out), self.ndays)
        for a, b in zip(self.traces, out):
            self.assertEqual(a.stats.starttime, b.stats.starttime)
            self.assertEqual(a.stats.npts, b.stats.npts)
        ref, _ = self.converter.apply(Stream([self.whole]), detrend=False)
        diff = (np.concatenate([t.data for t in out]).astype(int)
                - ref[0].data.astype(int))
        self.assertLessEqual(np.abs(diff).max(), 1)  # float rounding only
        # Inventory
        inv = read_inventory(str(self.root / 'stations_FIRconv.xml'))
        cha = inv.select(location='01', channel='LHZ')[0][0][0]
        stage = cha.response.response_stages[1]
        self.assertEqual(stage.symmetry, 'NONE')
        self.assertEqual(float(stage.decimation_delay), 0.)
        self.assertEqual(float(stage.decimation_correction), 0.)
        self.assertAlmostEqual(
            cha.response.instrument_sensitivity.value
            / inv.select(location='00')[0][0][0].response
            .instrument_sensitivity.value, 1.)

    def test_decimate_SDS(self):
        decimate_SDS.main(self._args(decimate_SDS, [
            str(self.root / 'SDS'), self.inv_file, str(self.sr), '5',
            '--out_dir', str(self.root / 'SDS_decim'), '-q']))
        out = self._files('SDS_decim/*/XX/STA/VHZ.D/XX.STA.00.*')
        self.assertEqual(len(out), self.ndays)
        decimator = Decimator([5])
        phase = decimate_SDS._grid_phase(self.traces[0], 5)
        ref = decimator.decimate(Stream([self.whole]), phase=phase)[0]
        self.assertEqual(out[0].stats.starttime, ref.stats.starttime)
        self.assertListEqual(
            np.concatenate([t.data for t in out]).tolist(), ref.data.tolist())
        inv = read_inventory(str(self.root / 'stations_decim.xml'))
        self.assertEqual(len(inv.select(channel='VHZ')), 1)
        # Re-running into the same output directory works
        decimate_SDS.main(self._args(decimate_SDS, [
            str(self.root / 'SDS'), self.inv_file, str(self.sr), '5',
            '--out_dir', str(self.root / 'SDS_decim'), '-q']))


def suite():
    return unittest.makeSuite(TestSDSScripts, 'test')


if __name__ == '__main__':
    unittest.main(defaultTest='suite')
