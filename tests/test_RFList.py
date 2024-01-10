#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Functions to test RFList (organises and names ResponseFunctions for DataCleaner)

Does RFList have any purpose, now that removed channel information is
stored in trace.stats?
"""
# from os import system
import unittest
import inspect
from pathlib import Path
import fnmatch

from tiskitpy import RFList, ResponseFunctions as RF

class TestMethods(unittest.TestCase):
    """
    Test suite
    """
    def setUp(self):
        self.path = Path(
            inspect.getfile(inspect.currentframe())).resolve().parent

    def _RFList_create(self):
        """
        Create a test RFList object

        Has no response function, just the channels information
        """
        channels = ['A.BH1', 'A.BH2', 'A.BHZ', 'A.BDH']
        rf_list = RFList()
        out_channels = channels.copy()
        for component in ("*BDH", '*BH1', '*BH2'):
            in_channel = fnmatch.filter(channels, component)[0]
            out_channels = [x for x in out_channels if not x == in_channel]
            rf = RF(None, in_channel, out_channels, quiet=True)
            rf_list.append(rf)
        return rf_list

    def test_RFList_channel_map(self):
        rf_list = self._RFList_create()
        self.assertEqual(
            f'{rf_list[0].__str__()}',
            'ResponseFunctions object:\n'
            "  input_channel_id='A.BDH'\n"
            "  output_channel_ids=['A.BH1', 'A.BH2', 'A.BHZ']\n"
            "  noise_channels=['output', 'output', 'output']\n"
            "  n_windows=0"
        )
        self.assertEqual(
            f'{rf_list[1].__str__()}',
            'ResponseFunctions object:\n'
            "  input_channel_id='A.BH1'\n"
            "  output_channel_ids=['A.BH2', 'A.BHZ']\n"
            "  noise_channels=['output', 'output']\n"
            "  n_windows=0"
        )
        self.assertEqual(
            f'{rf_list[2].__str__()}',
            'ResponseFunctions object:\n'
            "  input_channel_id='A.BH2'\n"
            "  output_channel_ids=['A.BHZ']\n"
            "  noise_channels=['output']\n"
            "  n_windows=0"
        )
        self.assertEqual(
            f'{rf_list.__str__()}',
            "RFList object:\n"
            "   Step | Channel to remove  | Channels to remove from\n"
            "   ==== | ================== | ========================================\n"
            "      1 | A.BDH              | ['A.BH1', 'A.BH2', 'A.BHZ']\n"
            "      2 | A.BH1              | ['A.BH2', 'A.BHZ']\n"
            "      3 | A.BH2              | ['A.BHZ']\n"
            )


def suite():
    return unittest.makeSuite(TestMethods, 'test')


if __name__ == '__main__':
    unittest.main(defaultTest='suite')
