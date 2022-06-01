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

from crawtools.spectral import (DCTFs, DCTF, remove_str, strip_remove_str,
                                TransferFunctions as TF)


class TestMethods(unittest.TestCase):
    """
    Test suite
    """
    def setUp(self):
        self.path = Path(
            inspect.getfile(inspect.currentframe())).resolve().parent

    def _DCTFs_create(self):
        """
        Create a test DCTF object

        Has no transfer function, just the channels information
        """
        channels = ['A.BH1', 'A.BH2', 'A.BHZ', 'A.BDH']
        dctfs = DCTFs()
        remove_seq = ''
        out_channels = channels.copy()
        for component in ("*BDH", '*BH1', '*BH2'):
            in_channel = fnmatch.filter(channels, component)[0]
            # print(f'{channels=}, {component=}, {channel=}')
            ic = in_channel + remove_seq
            remove_new = remove_str(in_channel, remove_seq)
            remove_seq += remove_new
            out_channels = [x + remove_seq
                            for x in strip_remove_str(out_channels)
                            if not x == strip_remove_str(in_channel)]
            tf = TF(None, ic, out_channels, quiet=True)
            dctfs.append(DCTF(ic, remove_seq, tf))
        return dctfs

    def test_DCTFs_channel_map(self):
        dctfs = self._DCTFs_create()
        mapper = dctfs._channel_map()
        self.assertTrue(mapper == {'A.BH1': 'A.BH1-H',
                                   'A.BH2': 'A.BH2-H-1',
                                   'A.BHZ': 'A.BHZ-H-1-2'})

    def test_DCTFs_update_channel_names(self):
        dctfs = self._DCTFs_create()
        channel_names = ['A.BH1', 'A.BH2', 'A.BHZ', 'A.BDH']
        new_channel_names = dctfs.update_channel_names(channel_names)
        self.assertTrue(new_channel_names == ['A.BH1-H', 'A.BH2-H-1',
                                              'A.BHZ-H-1-2', 'A.BDH'])
        channel_names = ['A.BH1', 'A.BH2', 'A.BDH'] # no BHZ
        self.assertRaises(ValueError, dctfs.update_channel_names,
                          channel_names)

    def test_DCTFs_update_channel_keys(self):
        dctfs = self._DCTFs_create()
        in_dict = {'A.BH1': 'hi', 'A.BH2': 'ho', 'A.BHZ': 'he', 'A.BDH': 'hu'}
        out_dict = dctfs.update_channel_keys(in_dict)
        # Verify new keys added
        self.assertTrue(out_dict['A.BH1-H'] == in_dict['A.BH1'])
        self.assertTrue(out_dict['A.BH2-H-1'] == in_dict['A.BH2'])
        self.assertTrue(out_dict['A.BHZ-H-1-2'] == in_dict['A.BHZ'])
        # Verify old keys kept
        self.assertTrue(out_dict['A.BH1'] == in_dict['A.BH1'])
        self.assertTrue(out_dict['A.BH2'] == in_dict['A.BH2'])
        self.assertTrue(out_dict['A.BHZ'] == in_dict['A.BHZ'])
        self.assertTrue(out_dict['A.BDH'] == in_dict['A.BDH'])


def suite():
    return unittest.makeSuite(TestMethods, 'test')


if __name__ == '__main__':
    unittest.main(defaultTest='suite')
