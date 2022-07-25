#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Functions to test spectral functions
"""
import unittest
# from pathlib import Path
# import inspect

from obspy.core.stream import read

from tiskit import SpectralDensity, DataCleaner
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
        self.dc = DataCleaner(self.stream, remove_list=["XX.STA.00.BX1"],
                              window_s=self.window_s)

    def test_str(self):
        """Test __str__ function"""
        # print(self.dc)
        self.assertEqual(
            self.dc.__str__(),
            "DataCleaner object:\n"
            "        Input channel | output channels\n"
            "   ================== | ===============\n"
            "   XX.STA.00.BX1      | ['XX.STA.00.BX2-1', 'XX.STA.00.BX3-1']\n",
        )

    def test_clean_sdf(self):
        """Test clean_sdf function"""
        sdf = SpectralDensity.from_stream(self.stream, window_s=self.window_s)
        for fast_calc in (False, True):
            cleaned = self.dc.clean_sdf(sdf, fast_calc=fast_calc)
            # cleaned.plot(overlay=True)
            # self.assertTrue(np.all(self.xf.freqs
            #                        == np.arange(1, num+1)*stop/num))

    def test_clean_stream_to_sdf(self):
        """Test clean_stream_to_sdf function"""
        for fast_calc in (False, True):
            # cleaned = self.dc.clean_stream_to_sdf(self.stream)
            cleaned = self.dc.clean_stream_to_sdf(self.stream,
                                                  window_s=self.window_s)
            # cleaned.plot(overlay=True)
            # self.assertTrue(np.all(self.xf.freqs
            #                 == np.arange(1, num+1)*stop/num))

    def test_clean_stream(self):
        """Test clean_stream function"""
        for itd in (False, True):
            cleaned = self.dc.clean_stream(self.stream, in_time_domain=itd)
            # cleaned.plot(end_time=cleaned[0].stats.starttime + 10)
            sdf_cleaned = SpectralDensity.from_stream(self.stream,
                                                      window_s=self.window_s)
        

def suite():
    return unittest.makeSuite(TestMethods, "test")


if __name__ == "__main__":
    unittest.main(defaultTest="suite")
