#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Functions to test spectral functions
"""
import unittest
from pathlib import Path
import inspect

from obspy.core.stream import Stream

from tiskitpy.utils import get_full_id, match_one_str
from make_test_stream import make_test_stream


class TestMethods(unittest.TestCase):
    """
    Test suite
    """

    def setUp(self):
        self.stream, _ = make_test_stream()

    def test_match_one_str(self):
        """Test add_trace()"""
        str_list = ['ABC', 'AFIT', 'DEF', 'GHI']
        self.assertEqual(match_one_str('ABC', str_list, 'str', 'list'), 'ABC')
        self.assertEqual(match_one_str('AB*', str_list, 'str', 'list'), 'ABC')
        self.assertEqual(match_one_str('AF*', str_list, 'str', 'list'), 'AFIT')
        self.assertEqual(match_one_str('*F', str_list, 'str', 'list'), 'DEF')
        self.assertEqual(match_one_str('??F', str_list, 'str', 'list'), 'DEF')
        with self.assertRaises(ValueError):
            match_one_str('A*', str_list, 'str', 'list')
        with self.assertRaises(ValueError):
            match_one_str('?F', str_list, 'str', 'list')
        with self.assertRaises(ValueError):
            match_one_str('B*', str_list, 'str', 'list')
        
    def test_get_full_id(self):
        """Test add_trace()"""
        self.assertEqual(get_full_id('*BX1', self.stream), 'XX.STA.00.BX1')
        self.assertEqual(get_full_id('*1', self.stream), 'XX.STA.00.BX1')
        self.assertEqual(get_full_id('*3', self.stream), 'XX.STA.00.BX3')
        self.assertEqual(get_full_id('*H', self.stream), 'XX.STA.00.BDH')
        with self.assertRaises(ValueError):
            get_full_id('*Z', self.stream)
        with self.assertRaises(ValueError):
            get_full_id('?1', self.stream)
        

def suite():
    return unittest.makeSuite(TestMethods, "test")


if __name__ == "__main__":
    unittest.main(defaultTest="suite")
