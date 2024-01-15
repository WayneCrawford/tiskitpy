#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Functions to test spectral functions
"""
# from os import system
import unittest
from pathlib import Path
from copy import deepcopy

from obspy.core import Stream, Trace

# sys.path.append("..")

from tiskitpy.utils import SeisRotate


class TestMethods(unittest.TestCase):
    """
    Test suite
    """
    def setUp(self):
        self.path = Path(__file__).resolve().parent
        self.test_path = self.path / "data" / "utils"

    def test_get_one_trace(self):
        tr = SeisRotate._get_one_trace(
            _quick_stream(['BH1', 'BH2', 'BH3', 'BHZ', 'BHP']), '1')
        self.assertEqual(tr.stats.channel, 'BH1')
        # Stream with missing component should raise IndexError
        with self.assertRaises(IndexError):
            SeisRotate._get_one_trace(_quick_stream(['BH1', 'BH2', 'BH3']), 'Z')
        # Stream with redundant component should raise ValueError
        with self.assertRaises(ValueError):
            SeisRotate._get_one_trace(_quick_stream(['BH1', 'BH2', 'SH2']), '2')

    def test_get_seis_traces(self):
        Z, N, E = SeisRotate._get_seis_traces(
            _quick_stream(['BH1', 'BH2', 'BH3', 'BHZ', 'BHP']))
        self.assertEqual(Z.stats.channel, 'BHZ')
        self.assertEqual(N.stats.channel, 'BH1')
        self.assertEqual(E.stats.channel, 'BH2')

    def test_separate_streams(self):
        seis, other = SeisRotate.separate_streams(
            _quick_stream(['BH1', 'BH2', 'BH3', 'BHZ', 'BHP']))
        self.assertEqual(len(seis), 3)
        self.assertEqual(len(other), 2)
        seis, other = SeisRotate.separate_streams(
            _quick_stream(['BHE', 'BHN', 'BHZ']))
        self.assertEqual(len(seis), 3)
        self.assertIsNone(other)
        # Stream with missing component should raise IndexError
        with self.assertRaises(IndexError):
            SeisRotate.separate_streams(
                _quick_stream(['BH1', 'BH2', 'BH3', 'BHP']))
        # Stream with redundant component should raise ValueError
        with self.assertRaises(ValueError):
            SeisRotate.separate_streams(_quick_stream(['BH1', 'BH2', 'BHZ', 'SH2']))

#     def test_cleaner_string_make(self):
#         self.assertEqual(CS.make('BHZ', '-X-Y'), '-Z')
#         self.assertEqual(CS.make('BHZ', '-X-Y', full_str=True), '-BHZ')
#         self.assertEqual(CS.make('BHX', '-X-Y'), '-HX')
#         with self.assertRaisesRegex(ValueError,
#                                     r"new_name\('X'\) already in prev_str\('-X-Y'\)"):
#             CS.make('X', '-X-Y')
# 
#     def test_cleaner_string_add(self):
#         self.assertEqual(CS.add('BHZ', '-X-Y'), '-X-Y-Z')
#         self.assertEqual(CS.add('BHX', '-X-Y'), '-X-Y-HX')
#         with self.assertRaisesRegex(ValueError,
#                                     r"new_name\('X'\) already in prev_str\('-X-Y'\)"):
#             CS.add('X', '-X-Y'),
#                                
#     def test_cleaner_string_insert(self):
#         self.assertEqual(CS.insert('-X-Y', 'BHZ'), '-X-Y.BHZ')
#         self.assertEqual(CS.insert('-X-Y', 'XX.STA.00.BHZ'), 'XX.STA.00-X-Y.BHZ')
#         self.assertEqual(CS.insert('-Z', 'XX.SSS.00.BDH'), 'XX.SSS.00-Z.BDH')
#         self.assertEqual(CS.insert('-ROT-X-Y', 'XX.SSS.00.BDH'),
#                          'XX.SSS.00-ROT-X-Y.BDH')
#         with self.assertRaisesRegex(ValueError,
#                                     r"cleaner_str \('-X'\) already in orig_id.location \('00-X-Y'\)"):
#             CS.insert('-X', 'XX.SSS.00-X-Y.BDH')
# 
#     def test_cleaner_string_insert_id(self):
#         self.assertEqual(CS.insert_id('BHX', 'BHZ'), '-X.BHZ')
#         self.assertEqual(CS.insert_id('BHX', 'BHZ', full_str=True), '-BHX.BHZ')
#         self.assertEqual(CS.insert_id('BHZ', 'XX.SSS.00.BDH'),
#                          'XX.SSS.00-Z.BDH')
#         self.assertEqual(CS.insert_id('BHZ', 'XX.SSS.00-X-Y.BDH'),
#                          'XX.SSS.00-X-Y-Z.BDH')
#         self.assertEqual(CS.insert_id('BHX', 'XX.SSS.00-X-Y.BDH'),
#                          'XX.SSS.00-X-Y-HX.BDH')
#         with self.assertRaisesRegex(ValueError,
#                                     r"new_name\('X'\) already in prev_str\('-X-Y'\)"):
#             CS.insert_id('X', 'XX.SSS.00-X-Y.BDH')
# 
#     def test_cleaner_string_strip(self):
#         self.assertEqual(CS.strip('-X-Y.BHZ'), 'BHZ')
#         self.assertEqual(CS.strip('.BHZ'), '.BHZ')
#         self.assertEqual(CS.strip('00-X-Y.BHZ'), '00.BHZ')
#         self.assertEqual(CS.strip('XX.STA.-X-Y.BHZ'), 'XX.STA..BHZ')
#         self.assertEqual(CS.strip('XX.STA.00.BHZ'), 'XX.STA.00.BHZ')
#         self.assertEqual(CS.strip('XX.SSS.00.BHZ'), 'XX.SSS.00.BHZ')
#         self.assertEqual(CS.strip('XX.SSS.00-X-Y.BHZ'), 'XX.SSS.00.BHZ')
#         self.assertEqual(CS.strip(['XX.SSS.00-1-2-H.BHZ', 'XX.SSS.00-1-2.BDH',
#                                    'XX.SSS.00-1.BH2', 'XX.SSS.00.BH1']),
#                          ['XX.SSS.00.BHZ', 'XX.SSS.00.BDH',
#                           'XX.SSS.00.BH2', 'XX.SSS.00.BH1'])
# 
#     def test_cleaner_string_extract(self):
#         self.assertEqual(CS.extract('XX.SSS.00.BHZ'), '')
#         self.assertEqual(CS.extract('XX.SSS.00-X-Y.BHZ'), '-X-Y')
#         self.assertEqual(CS.extract('00-X-Y.BHZ'), '-X-Y')
#         self.assertEqual(CS.extract(['XX.SSS.00-1-2-H.BHZ', 'XX.SSS.00-1-2.BDH',
#                                      'XX.SSS.00-1.BH2', 'XX.SSS.00.BH1']),
#                          ['-1-2-H', '-1-2', '-1', ''])
# 
#     def test_cleaner_string_separate(self):
#         self.assertEqual(CS._separate('XX.SSS.00.BHZ'), ['XX.SSS.00.BHZ', ''])
#         self.assertEqual(CS._separate('XX.SSS.00-X-Y.BHZ'), ['XX.SSS.00.BHZ', '-X-Y'])
#         self.assertEqual(CS._separate('00-X-Y.BHZ'), ['00.BHZ', '-X-Y'])
#         tester = ['XX.SSS.00-1-2-H.BHZ', 'XX.SSS.00-1-2.BDH',
#                   'XX.SSS.00-1.BH2', 'XX.SSS.00.BH1']
#         tester_copy = deepcopy(tester)
#         self.assertEqual(CS._separate(tester),
#                          [['XX.SSS.00.BHZ', 'XX.SSS.00.BDH',
#                            'XX.SSS.00.BH2', 'XX.SSS.00.BH1'],
#                            ['-1-2-H', '-1-2', '-1', '']])
#         self.assertEqual(tester, tester_copy)
# 
#     def test_remove_cleaner_string(self):
#         orig_id = 'XX.SSS.00.BHZ'
#         cleaner_string = '-ROT-X-Y'
#         clean_id = CS.insert(cleaner_string, orig_id)
#         self.assertEqual(clean_id, 'XX.SSS.00-ROT-X-Y.BHZ')
#         # test string input
#         self.assertEqual(remove_cleaner_string(clean_id), orig_id)
#         # test Stream input
#         stream = _quick_stream(('BHX', 'BHY', 'BHZ'))
#         for tr in stream:
#             tr.stats.network = 'XX'
#             tr.stats.station = 'SSS'
#             tr.stats.location = '00'
#             tr.id = CS.insert_id('BDH', tr.get_id())
#             self.assertEqual(tr.get_id(), f'XX.SSS.00-H.{tr.stats.channel}')
#         new_stream = remove_cleaner_string(stream)
#         for tr in new_stream:
#             self.assertEqual(tr.get_id(), f'XX.SSS.00.{tr.stats.channel}')
#         # test Trace input
#         tr = _quick_stream(('BHX',))[0]
#         tr.stats.network = 'XX'
#         tr.stats.station = 'SSS'
#         tr.stats.location = '00'
#         tr.id = CS.insert_id('BDH', tr.get_id())
#         self.assertEqual(tr.get_id(), f'XX.SSS.00-H.BHX')
#         tr = remove_cleaner_string(tr)
#         self.assertEqual(tr.get_id(), f'XX.SSS.00.BHX')
        
def _quick_stream(chan_list):
    "Stream with given channels"
    stream = Stream()
    for chan in chan_list:
        stream += Trace(header={"channel": chan})
    return stream


def suite():
    return unittest.makeSuite(TestMethods, 'test')


if __name__ == '__main__':
    unittest.main(defaultTest='suite')
