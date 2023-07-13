#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Functions to test spectral functions
"""
# from os import system
import unittest
import inspect
from pathlib import Path
import fnmatch

from tiskitpy import DCRFs, DCRF, CleanerString as CS, ResponseFunctions as RF


class TestMethods(unittest.TestCase):
    """
    Test suite
    """
    def setUp(self):
        self.path = Path(
            inspect.getfile(inspect.currentframe())).resolve().parent

    def _DCRFs_create(self):
        """
        Create a test DCRF object

        Has no response function, just the channels information
        """
        channels = ['A.BH1', 'A.BH2', 'A.BHZ', 'A.BDH']
        dcrfs = DCRFs()
        remove_seq = ''
        out_channels = channels.copy()
        for component in ("*BDH", '*BH1', '*BH2'):
            in_channel = fnmatch.filter(channels, component)[0]
            # print(f'{channels=}, {component=}, {channel=}')
            ic = in_channel + remove_seq
            remove_new = CS.make(in_channel, remove_seq)
            remove_seq += remove_new
            # out_channels = [x + remove_seq
            out_channels = [CS.insert(remove_seq, x)
                            for x in CS.strip(out_channels)
                            if not x == CS.strip(in_channel)]
            rf = RF(None, ic, out_channels, quiet=True)
            dcrfs.append(DCRF(ic, remove_seq, rf))
        return dcrfs

    def test_DCRFs_channel_map(self):
        dcrfs = self._DCRFs_create()
        mapper = dcrfs._channel_map()
        self.assertTrue(mapper == {'A.BH1': 'A-H.BH1',
                                   'A.BH2': 'A-H-1.BH2',
                                   'A.BHZ': 'A-H-1-2.BHZ'})

    def test_DCRFs_update_channel_names(self):
        dcrfs = self._DCRFs_create()
        channel_names = ['A.BH1', 'A.BH2', 'A.BHZ', 'A.BDH']
        new_channel_names = dcrfs.update_channel_names(channel_names)
        # self.assertTrue(new_channel_names == ['A.BH1-H', 'A.BH2-H-1',
        #                                       'A.BHZ-H-1-2', 'A.BDH'])
        self.assertTrue(new_channel_names == ['A-H.BH1', 'A-H-1.BH2',
                                              'A-H-1-2.BHZ', 'A.BDH'])
        channel_names = ['A.BH1', 'A.BH2', 'A.BDH']  # no BHZ
        self.assertRaises(ValueError, dcrfs.update_channel_names,
                          channel_names)

    def test_DCRFs_update_channel_keys(self):
        dcrfs = self._DCRFs_create()
        in_dict = {'A.BH1': 'hi', 'A.BH2': 'ho', 'A.BHZ': 'he', 'A.BDH': 'hu'}
        out_dict = dcrfs.update_channel_keys(in_dict)
        # Verify new keys added
        self.assertTrue(out_dict['A-H.BH1'] == in_dict['A.BH1'])
        self.assertTrue(out_dict['A-H-1.BH2'] == in_dict['A.BH2'])
        self.assertTrue(out_dict['A-H-1-2.BHZ'] == in_dict['A.BHZ'])
        # Verify old keys kept
        self.assertTrue(out_dict['A.BH1'] == in_dict['A.BH1'])
        self.assertTrue(out_dict['A.BH2'] == in_dict['A.BH2'])
        self.assertTrue(out_dict['A.BHZ'] == in_dict['A.BHZ'])
        self.assertTrue(out_dict['A.BDH'] == in_dict['A.BDH'])


def suite():
    return unittest.makeSuite(TestMethods, 'test')


if __name__ == '__main__':
    unittest.main(defaultTest='suite')