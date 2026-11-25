#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Test the SDS console scripts (decimate_SDS, FIR_convert_SDS, compare_SDS_XML)
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
from tiskitpy.console_scripts import decimate_SDS, fir_convert, compare_SDS_XML


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

    def _add_sds_channel(self, chan, loc='00'):
        """Add a one-file channel to the SDS, return its path"""
        tr = self.traces[0].copy()
        tr.stats.channel = chan
        tr.stats.location = loc
        fpath = (self.root / 'SDS' / '2015' / 'XX' / 'STA' / f'{chan}.D'
                 / f'XX.STA.{loc}.{chan}.D.2015.365')
        fpath.parent.mkdir(parents=True, exist_ok=True)
        tr.write(str(fpath), format='MSEED')
        return fpath

    def test_FIR_convert_SDS_copy(self):
        # A channel without the FIR and a channel not in the StationXML,
        # which should both be copied to out_dir
        no_fir = self._add_sds_channel('LHN')
        sds_only = self._add_sds_channel('LHE')
        inv = read_inventory(self.inv_file)
        extra = inv[0][0].channels[0].copy()
        extra.code = 'LHN'
        extra.response.response_stages = extra.response.response_stages[:1]
        inv[0][0].channels.append(extra)
        inv.write(self.inv_file, format='STATIONXML')
        out_dir = self.root / 'SDS_conv'
        base = [str(self.root / 'SDS'), self.inv_file,
                '--builtin', 'lc2000_fir3_0', '-q', '--out_dir']
        fir_convert.main(self._args(fir_convert, base + [str(out_dir)]))
        self.assertEqual(
            len(list(out_dir.glob('*/XX/STA/LHZ.D/XX.STA.01.*'))), self.ndays)
        for fpath in (no_fir, sds_only):
            copied = out_dir / fpath.relative_to(self.root / 'SDS')
            self.assertTrue(copied.is_file())
            self.assertEqual(copied.read_bytes(), fpath.read_bytes())
        # The converted channel's source is not copied
        self.assertEqual(list(out_dir.glob('*/XX/STA/LHZ.D/XX.STA.00.*')),
                         [])
        # out_dir must not exist, unless --out_dir_exists_ok
        with self.assertRaises(ValueError):
            fir_convert.main(self._args(fir_convert, base + [str(out_dir)]))
        # With --out_dir_exists_ok, existing NSLCs are not overwritten
        copied = out_dir / no_fir.relative_to(self.root / 'SDS')
        copied.write_bytes(b'not overwritten')
        with self.assertLogs('tiskitpy', 'INFO') as cm:
            fir_convert.main(self._args(fir_convert, base + [
                str(out_dir), '--out_dir_exists_ok', '--loc_code', '02']))
        self.assertEqual(
            len(list(out_dir.glob('*/XX/STA/LHZ.D/XX.STA.02.*'))), self.ndays)
        self.assertEqual(copied.read_bytes(), b'not overwritten')
        summary = self._summary(cm)
        self.assertEqual(summary['XX.STA.00.LHZ'], 'converted')
        self.assertEqual(summary['XX.STA.00.LHN'], 'aborted')
        self.assertEqual(summary['XX.STA.00.LHE'], 'aborted')
        # Same location code: the conversion is aborted
        with self.assertLogs('tiskitpy', 'INFO') as cm:
            fir_convert.main(self._args(fir_convert, base + [
                str(out_dir), '--out_dir_exists_ok']))
        self.assertEqual(self._summary(cm)['XX.STA.00.LHZ'], 'aborted')
        # --no_copy
        out_dir2 = self.root / 'SDS_conv2'
        fir_convert.main(self._args(fir_convert, base + [str(out_dir2),
                                                         '--no_copy']))
        self.assertEqual(list(out_dir2.glob('*/XX/STA/LH[NE].D/*')), [])

    @staticmethod
    def _summary(cm):
        """{seed_id: action} from the summary in captured log output"""
        msg = [r.getMessage() for r in cm.records
               if r.getMessage().startswith('Summary: ')][-1]
        return {line.split(':')[0].strip(): line.split(':')[1].split()[0]
                for line in msg.splitlines()[1:]}

    def test_FIR_convert_SDS_existing_nslc(self):
        # The converted channel's NSLC already exists in the SDS
        existing = self._add_sds_channel('LHZ', loc='01')
        before = existing.read_bytes()
        # In place: the conversion is aborted
        with self.assertLogs('tiskitpy', 'INFO') as cm:
            inv = fir_convert.FIR_convert_SDS(self._args(fir_convert, [
                str(self.root / 'SDS'), self.inv_file,
                '--builtin', 'lc2000_fir3_0', '-q']),
                read_inventory(self.inv_file))
        summary = self._summary(cm)
        self.assertEqual(summary['XX.STA.00.LHZ'], 'aborted')
        self.assertEqual(summary['XX.STA.01.LHZ'], 'skipped')
        self.assertEqual(existing.read_bytes(), before)
        self.assertEqual(len(inv.select(location='01')), 0)
        # To out_dir: XX.STA.01.LHZ is copied, so the conversion is aborted
        out_dir = self.root / 'SDS_conv'
        with self.assertLogs('tiskitpy', 'INFO') as cm:
            fir_convert.main(self._args(fir_convert, [
                str(self.root / 'SDS'), self.inv_file, '--builtin',
                'lc2000_fir3_0', '-q', '--out_dir', str(out_dir)]))
        summary = self._summary(cm)
        self.assertEqual(summary['XX.STA.00.LHZ'], 'aborted')
        # ...and as nothing is converted, nothing is copied
        self.assertEqual(summary['XX.STA.01.LHZ'], 'skipped')
        self.assertFalse(out_dir.exists())

    def test_FIR_convert_SDS_nothing_to_convert(self):
        # No channel has the FIR: nothing is converted or copied
        inv = read_inventory(self.inv_file)
        resp = inv[0][0].channels[0].response
        resp.response_stages = resp.response_stages[:1]
        inv.write(self.inv_file, format='STATIONXML')
        out_dir = self.root / 'SDS_conv'
        fir_convert.main(self._args(fir_convert, [
            str(self.root / 'SDS'), self.inv_file, '--builtin',
            'lc2000_fir3_0', '-q', '--out_dir', str(out_dir)]))
        self.assertFalse(out_dir.exists())

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

    def test_compare_SDS_XML(self):
        # Add a channel only in the SDS and one only in the StationXML
        tr = self.traces[0].copy()
        tr.stats.channel = 'LHN'
        fpath = (self.root / 'SDS' / '2015' / 'XX' / 'STA' / 'LHN.D'
                 / 'XX.STA.00.LHN.D.2015.365')
        fpath.parent.mkdir(parents=True, exist_ok=True)
        tr.write(str(fpath), format='MSEED')
        inv = read_inventory(self.inv_file)
        sta = inv[0][0]
        extra = sta.channels[0].copy()
        extra.code = 'LHE'
        sta.channels.append(extra)
        sta.channels.append(sta.channels[0].copy())  # second epoch
        inv.write(self.inv_file, format='STATIONXML')
        comp = compare_SDS_XML.compare_SDS_XML(str(self.root / 'SDS'),
                                               self.inv_file)
        self.assertEqual(comp['both'], [('XX', 'STA', '00', 'LHZ')])
        self.assertEqual(comp['inv_only'], [('XX', 'STA', '00', 'LHE')])
        self.assertEqual(comp['sds_only'], [('XX', 'STA', '00', 'LHN')])


def suite():
    return unittest.makeSuite(TestSDSScripts, 'test')


if __name__ == '__main__':
    unittest.main(defaultTest='suite')
