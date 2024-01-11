#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Functions to test spectral functions
"""
import unittest

from obspy.core.stream import read

from tiskitpy import SpectralDensity, DataCleaner
from make_test_stream import make_test_stream


class TestMethods(unittest.TestCase):
    """
    Test suite
    """

    def setUp(self):
        # self.path = (Path(inspect.getfile(inspect.currentframe()))
        #              .resolve().parent)
        # self.test_path = self.path / "data" / "data_cleaner"
        self.window_s = 100.0

        self.stream, self.sineparms = make_test_stream()
        self.dc1 = DataCleaner(self.stream, remove_list=["XX.STA.00.BX1"],
                               window_s=self.window_s)
        self.dc12 = DataCleaner(self.stream, remove_list=["*1", "*2"],
                                window_s=self.window_s)

    def test_str(self):
        """Test __str__ function"""
        self.assertEqual(
            self.dc1.__str__(),
            "DataCleaner object:\n"
            " 6 starttimes\n"
            " response_functions = RFList object:\n"
            "   Step | Channel to remove  | Channels to remove from\n"
            "   ==== | ================== | ========================================\n"
            "      1 | XX.STA.00.BX1      | ['XX.STA.00.BX2', 'XX.STA.00.BX3', 'XX.STA.00.BDH']\n"
        )
        self.maxDiff=None
        self.assertEqual(
            self.dc12.__str__(),
            "DataCleaner object:\n"
            " 6 starttimes\n"
            " response_functions = RFList object:\n"
            "   Step | Channel to remove  | Channels to remove from\n"
            "   ==== | ================== | ========================================\n"
            "      1 | XX.STA.00.BX1      | ['XX.STA.00.BX2', 'XX.STA.00.BX3', 'XX.STA.00.BDH']\n"
            "      2 | XX.STA.00-1.BX2    | ['XX.STA.00-1.BX3', 'XX.STA.00-1.BDH']\n"
        )

    def test_clean_sdf(self):
        """Test clean_sdf function"""
        sdf = SpectralDensity.from_stream(self.stream, window_s=self.window_s)
        seed_ids = ['XX.STA.00.BX1', 'XX.STA.00.BX2', 'XX.STA.00.BX3', 'XX.STA.00.BDH']
        ch_ids = ['XX.STA.00.BX1', 'XX.STA.00-1.BX2', 'XX.STA.00-1.BX3', 'XX.STA.00-1.BDH']
        
        # Test dc1 (one channel removed)
        cleaned = self.dc1.clean_sdf(sdf)
        self.assertEqual(cleaned.ids, ch_ids)
        self.assertEqual(cleaned.seed_ids, seed_ids)
        for ch_id, cleaned_ids in {
            'XX.STA.00.BX1': None,
            'XX.STA.00-1.BX2': ['XX.STA.00.BX1'],
            'XX.STA.00-1.BX3': ['XX.STA.00.BX1'],
            'XX.STA.00-1.BDH': ['XX.STA.00.BX1']}.items():
            self.assertEqual(cleaned.clean_sequence(chan), cleaned_chans)
        
        # Test dc1 (two channels removed)
        ch_ids = ['XX.STA.00.BX1', 'XX.STA.00-1.BX2', 'XX.STA.00-1-2.BX3', 'XX.STA.00-1-2.BDH']
        cleaned = self.dc12.clean_sdf(sdf)
        self.assertEqual(cleaned.ids, ch_ids)
        self.assertEqual(cleaned.seed_ids, seed_ids)
        for ch_id, cleaned_ids in {
            'XX.STA.00.BX1': None,
            'XX.STA.00-1.BX2': ['XX.STA.00.BX1'],
            'XX.STA.00-1-2.BX3': ['XX.STA.00.BX1','XX.STA.00-1.BX2'],
            'XX.STA.00-1-2.BDH': ['XX.STA.00.BX1','XX.STA.00-1.BX2']}.items():
            self.assertEqual(cleaned.clean_sequence(chan), cleaned_chans)

    def test_clean_stream_to_sdf(self):
        """Test clean_stream_to_sdf function"""
        for fast_calc in (False, True):
            # cleaned = self.dc.clean_stream_to_sdf(self.stream)
            cleaned = self.dc1.clean_stream_to_sdf(self.stream,
                                                   window_s=self.window_s)
            # cleaned.plot(overlay=True)

    def test_clean_stream(self):
        """Test clean_stream function"""
        # for itd in (False, True):
        # Time domain too slow for these tests
        clean_tests = [{'cleaner': self.dc1,
                        'relations': {'XX.STA.00.BX1': None,
                                      'XX.STA.00.BX2': ['XX.STA.00.BX1'],
                                      'XX.STA.00.BX3': ['XX.STA.00.BX1'],
                                      'XX.STA.00.BDH': ['XX.STA.00.BX1'],
                                     }
                       },
                       {'cleaner': self.dc12,
                        'relations': {'XX.STA.00.BX1': None,
                                      'XX.STA.00.BX2': ['XX.STA.00.BX1'],
                                      'XX.STA.00.BX3': ['XX.STA.00.BX1', 'XX.STA.00.BX2'],
                                      'XX.STA.00.BDH': ['XX.STA.00.BX1', 'XX.STA.00.BX2'],
                                     }
                       }
                   ]
        for itd in (False,):
            for dc in clean_tests:
                cleaner = dc['cleaner']
                relations = dc['relations']
                cleaned = cleaner.clean_stream(self.stream, in_time_domain=itd)
                self.assertEqual([x.id for x in cleaned], list(relations.keys()))
                for x in cleaned:
                    self.assertEqual(x.stats.get('clean_sequence', None),
                                     relations[x.id])
                # cleaned.plot(end_time=cleaned[0].stats.starttime + 10)
                sdf_cleaned = SpectralDensity.from_stream(
                    self.stream, data_cleaner=cleaner, window_s=self.window_s)
                self.assertEqual(sdf_cleaned.ids, list(relations.keys()))
                for chan, clean_seq in relations.items():
                    self.assertEqual(sdf_cleaned.clean_sequence(chan),
                                     clean_seq)

def suite():
    return unittest.makeSuite(TestMethods, "test")


if __name__ == "__main__":
    unittest.main(defaultTest="suite")
