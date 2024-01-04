#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Functions to test spectral functions
"""
import unittest
from pathlib import Path
import inspect

from obspy.core.stream import read

from tiskitpy import CleanSequence as CS
from make_test_stream import make_test_stream


class TestMethods(unittest.TestCase):
    """
    Test suite
    """

    def setUp(self):
        self.path = (Path(inspect.getfile(inspect.currentframe()))
                     .resolve().parent)
        self.test_path = self.path / "data" / "clean_sequence"
        self.stream, _ = make_test_stream()
        self.trace = self.stream.select(channel='*3')[0].copy()

    def test_min_codes(self):
        """
        Test the private function behind _prepend_unique_str()
        """
        self.assertEqual(CS._min_codes(['BDH', 'BH1', 'BH2']),
                         ['H', '1', '2'])
        self.assertEqual(CS._min_codes(['BDH', 'BH1', 'BDH']),
                         ['H', '1', 'H'])
        self.assertEqual(CS._min_codes(['BDH', 'BH1', 'BCH']),
                         ['DH', '1', 'CH'])
        self.assertEqual(CS._min_codes(['01', '00', '00']),
                         ['1', '0', '0'])

    def test_prepend_unique_str(self):
        """
        Test the private function behind _id_list_str()
        """
        self.assertEqual(CS._prepend_unique_str(['BDH', 'BH1', 'BH2'],
                                                [None, None, None],
                                                [False, False, False]),
                         (['H', '1', '2'], [True, True, True]))
        self.assertEqual(CS._prepend_unique_str(['BDH', 'BH1', 'BH2'],
                                                [None, None, None],
                                                [False, False, False],
                                                'min_code'),
                         (['BDH', 'BH1', 'BH2'], [True, True, True]))
        self.assertEqual(CS._prepend_unique_str(['BDH', 'BH1', 'BH2'],
                                                ["", "", ""],
                                                [False, False, False]),
                         (['H_', '1_', '2_'], [True, True, True]))
        self.assertEqual(CS._prepend_unique_str(['BDH', 'BH1', 'BDH'],
                                                [None, None, None],
                                                [False, False, False]),
                         (['H', '1', 'H'], [False, True, False]))
        self.assertEqual(CS._prepend_unique_str(['01', '00', '00'],
                                                ["H", "1", "H"],
                                                [False, True, False]),
                         (['1_H', '1', '0_H'], [True, True, True]))
        self.assertEqual(CS._prepend_unique_str(['01', '00', '00'],
                                                ["H", "1", "H"],
                                                [False, True, False],
                                                'min_code'),
                         (['01_H', '1', '00_H'], [True, True, True]))

    def test_id_list_str(self):
        """
        Test the private function that creates a cleaner string
        """
        self.assertEqual(CS._id_list_str(['NN.SSS.LL.BDH',
                                          'NN.SSS.LL.BH1',
                                          'NN.SSS.LL.BH2']),
                         '-H-1-2')
        self.assertEqual(CS._id_list_str(['NN.SSS.LL.BDH',
                                          'NN.SSS.LL.BH1',
                                          'NN.SSS.LL.BH2'],
                                          'min_code'),
                         '-BDH-BH1-BH2')
        self.assertEqual(CS._id_list_str(['NN.SSS.LL.BDH',
                                          'NN.SSS.LL.BH1',
                                          'NN.SSS.LL.BH2'],
                                          'full'),
                         '-NN_SSS_LL_BDH-NN_SSS_LL_BH1-NN_SSS_LL_BH2')
        # Same channel removed twice raises ValueError
        with self.assertRaises(ValueError):
            CS._id_list_str(['NN.SSS.LL.BDH',
                             'NN.SSS.LL.BH1',
                             'NN.SSS.LL.BDH',
                             'NN.SSS.LL.BH2'])
        self.assertEqual(CS._id_list_str(['NN.SSS.L1.BDH',
                                          'NN.SSS.L1.BH1',
                                          'NN.SSS.L2.BDH',
                                          'NN.SSS.L1.BH2']),
                         '-1_H-1-2_H-2')
        self.assertEqual(CS._id_list_str(['NN.SSS.L1.BDH',
                                          'NN.SSS.L1.BH1',
                                          'NN.SSS.L2.BDH',
                                          'NN.SSS.L1.BH2'],
                                          'min_code'),
                         '-L1_BDH-BH1-L2_BDH-BH2')
        self.assertEqual(CS._id_list_str(['NN.SSS.L1.BDH',
                                          'NN.SSS.L1.BH1',
                                          'NN.SSS.L2.BDH',
                                          'NN.SSS.L1.BH2'],
                                          'min_level'),
                         '-L1_BDH-L1_BH1-L2_BDH-L1_BH2')

    def test_tag(self):
        cleaned_trace = CS.tag(self.trace, '...BDH')
        cleaned_stream = self.stream.copy()
        cleaned_stream = CS.tag(cleaned_stream, '...BDH', cleaned_ids=('*BX1', '*BX2', '*BX3'))
        cleaned_stream = CS.tag(cleaned_stream, '...BX1', cleaned_ids=('*BX2', '*BX3'))
        cleaned_stream = CS.tag(cleaned_stream, '...BX2', cleaned_ids=('*BX3',))
        # Test trace tag()
        new_trace = CS.tag(self.trace.copy(), 'BDH')
        self.assertEqual(new_trace.get_id(), 'XX.STA.00.BX3')
        self.assertEqual(CS.seedid_tag(new_trace).get_id(), 'XX.STA.00-H.BX3')
        self.assertEqual(CS.seedid_untag(CS.seedid_tag(new_trace)).get_id(), 'XX.STA.00.BX3')
        # Test trace tag() with a list of cleaned channels
        new_trace = CS.tag(self.trace.copy(), ['RRR', 'BH1', 'BDH'])
        self.assertEqual(new_trace.get_id(), 'XX.STA.00.BX3')
        self.assertEqual(CS.seedid_tag(new_trace).get_id(), 'XX.STA.00-R-1-H.BX3')
        self.assertEqual(CS.seedid_untag(CS.seedid_tag(new_trace)).get_id(), 'XX.STA.00.BX3')
        
        # Now test stream tag()
        test_stream = self.stream.copy()
        # Add same cleaned channel tag to all traces in a stream
        new_test_stream = CS.tag(test_stream, 'BX3')
        for c in ('BDH', 'BX1', 'BX2', 'BX3'):
            x = new_test_stream.select(channel=c)[0]
            self.assertEqual(x.get_id(), f'XX.STA.00.{c}')
            self.assertEqual(CS.seedid_tag(x).get_id(), f'XX.STA.00-3.{c}')

        # Add a different cleaner str to each trace, using stream interface
        test_stream = self.stream.copy()
        test_stream = CS.tag(test_stream, 'BDH', ('BX1', 'BX2', 'BX3'))
        test_stream = CS.tag(test_stream, 'BX1', ('BX2', 'BX3'))
        test_stream = CS.tag(test_stream, 'BX2', ('BX3',))
        self.assertEqual(test_stream.select(channel='BDH')[0].get_id(), 'XX.STA.00.BDH')
        self.assertEqual(test_stream.select(channel='BX1')[0].get_id(), 'XX.STA.00.BX1')
        self.assertEqual(test_stream.select(channel='BX2')[0].get_id(), 'XX.STA.00.BX2')
        self.assertEqual(test_stream.select(channel='BX3')[0].get_id(), 'XX.STA.00.BX3')
        self.assertEqual(CS.seedid_tag(test_stream.select(channel='BX1')[0]).get_id(), 'XX.STA.00-H.BX1')
        self.assertEqual(CS.seedid_tag(test_stream.select(channel='BX2')[0]).get_id(), 'XX.STA.00-H-1.BX2')
        self.assertEqual(CS.seedid_tag(test_stream.select(channel='BX3')[0]).get_id(), 'XX.STA.00-H-1-2.BX3')
        
        # Save the cleaned trace and stream for further tests
        self.cleaned_trace = cleaned_trace
        self.cleaned_stream = cleaned_stream
# 
#     def test_remove_and_add_loc(self):
#         """Test remove_loc() and put_loc()"""
#         test_trace = self.cleaned_trace.copy()
#         self.assertEqual(test_trace.get_id(), 'XX.STA.00-BDH.BX3')
#         self.assertEqual(CS.del_id(test_trace).get_id(), 'XX.SS.00.BX3')
#         self.assertEqual(CS.put_id(test_trace).get_id(), 'XX.SS.00-BDH.BX3')
#         
#         # Now test stream operations
#         self.assertEqual(self.cleaned_stream.select(channel='BDH').get_id(), 'XX.STA.00.BDH')
#         self.assertEqual(self.cleaned_stream.select(channel='BX1').get_id(), 'XX.STA.00-H.BX1')
#         self.assertEqual(self.cleaned_stream.select(channel='BX2').get_id(), 'XX.STA.00-H-1.BX2')
#         self.assertEqual(self.cleaned_stream.select(channel='BX3').get_id(), 'XX.STA.00-H-1-2.BX3')
#         
#         # Remove clean_string from location code
#         unloc_stream = CS.del_id(self.cleaned_stream)
#         self.assertEqual(unloc_stream.select(channel='BDH').get_id(), 'XX.STA.00.BDH')
#         self.assertEqual(unloc_stream.select(channel='BX1').get_id(), 'XX.STA.00.BX1')
#         self.assertEqual(unloc_stream.select(channel='BX2').get_id(), 'XX.STA.00.BX2')
#         self.assertEqual(unloc_stream.select(channel='BX3').get_id(), 'XX.STA.00.BX3')
# 
#         # Put clean_strings back in
#         loc_stream = CS.put_id(self.cleaned_stream)
#         self.assertEqual(unloc_stream.select(channel='BDH').get_id(), 'XX.STA.00.BDH')
#         self.assertEqual(unloc_stream.select(channel='BX1').get_id(), 'XX.STA.00-H.BX1')
#         self.assertEqual(unloc_stream.select(channel='BX2').get_id(), 'XX.STA.00-H-1.BX2')
#         self.assertEqual(unloc_stream.select(channel='BX3').get_id(), 'XX.STA.00-H-1-2.BX3')
# 
#         # doesn't double the loc stream, gives warning
#         loc_stream = CS.put_id(self.loc_stream)
#         self.assertEqual(unloc_stream.select(channel='BDH').get_id(), 'XX.STA.00.BDH')
#         self.assertEqual(unloc_stream.select(channel='BX1').get_id(), 'XX.STA.00-H.BX1')
#         self.assertEqual(unloc_stream.select(channel='BX2').get_id(), 'XX.STA.00-H-1.BX2')
#         self.assertEqual(unloc_stream.select(channel='BX3').get_id(), 'XX.STA.00-H-1-2.BX3')
                
#     def test_string(self):
#         """Test string() (returns cleaner string)"""
#         self.assertEqual(CS.string(self.cleaned_trace), '-BDH')
#         self.assertEqual(self.cleaned_trace.stats.clean_sequence, ['BDH'])
#         self.assertEqual(CS.string(self.cleaned_stream), ['', '-H', '-H-1', '-H-1-2'])
# 
#     def test_str_make(self):
#         """Test make()"""
#         self.assertEqual(CS.str_make('BHZ', '-X-Y'), '-Z')
#         self.assertEqual(CS.str_make('BHZ', '-X-Y', full_str=True), '-BHZ')
#         self.assertEqual(CS.str_make('BHX', '-X-Y'), '-HX')
#         with self.assertRaises(ValueError):
#             CS.str_make('X', '-X-Y')
# 
#     def test_str_insert(self):
#         """Test insert()"""
#         self.assertEqual(CS.str_insert('-X-Y', 'BHZ'), '-X-Y.BX3')
#         self.assertEqual(CS.str_insert('-X-Y', 'XX.STA.00.BHZ'), 'XX.STA.00-X-Y.BHZ')
#         
#     def test_str_insert_id(self):
#         """Test insert_id()"""
#         self.assertEqual(CS.str_insert_id('BHX', 'BHZ'), '-X.BHZ')
#         self.assertEqual(CS.str_insert_id('BHX', 'BHZ', full_str=True), '-BHX.BHZ')
#         self.assertEqual(CS.str_insert_id('BHX', 'XX.STA.00.BHZ'), 'XX.STA.00-X.BHZ')
#         self.assertEqual(CS.str_insert_id('BHX', 'XX.STA.00-Y.BHZ'), 'XX.STA.00-Y-X.BHZ')
#         
#     def test_str_strip(self):
#         """Test strip()"""
#         self.assertEqual(CS.str_strip('-X-Y.BHZ'), 'BHZ')
#         self.assertEqual(CS.str_strip('00-X-Y.BHZ'), '00.BHZ')
#         self.assertEqual(CS.str_strip('.BHZ'), '.BHZ')
#         self.assertEqual(CS.str_strip('XX.STA.00-X-Y.BHZ'), 'XX.STA.00.BHZ')
#         self.assertEqual(CS.str_strip('XX.STA.00-X-Y.BHZ'), 'XX.STA.00.BHZ')
#         self.assertEqual(CS.str_strip('XX.STA.00.BHZ'), 'XX.STA.00.BHZ')
        

def suite():
    return unittest.makeSuite(TestMethods, "test")


if __name__ == "__main__":
    unittest.main(defaultTest="suite")
