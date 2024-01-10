#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Functions to test CleanedStream class
"""
import unittest

from obspy.core.stream import Stream

from tiskitpy import CleanedStream, CleanSequence as CS
from make_test_stream import make_test_stream


class TestMethods(unittest.TestCase):
    """
    Test suite
    """

    def setUp(self):
        self.stream, _ = make_test_stream()
        self.stream = CS.tag(self.stream, '...BDH', cleaned_ids=('*BX1', '*BX2', '*BX3'))
        self.stream = CS.tag(self.stream, '...BX1', cleaned_ids=('*BX2', '*BX3'))
        self.stream = CS.tag(self.stream, '...BX2', cleaned_ids=('*BX3',))
        self.maxDiff = None

    def test_str(self):
        """Test __str__ function"""
        str_standard = (
            "4 Trace(s) in Stream:\n"
            "XX.STA.00.BDH | 2020-01-01T00:00:00.000000Z - 2020-01-01T00:16:39.990000Z | 100.0 Hz, 100000 samples\n"
            "XX.STA.00.BX1 | 2020-01-01T00:00:00.000000Z - 2020-01-01T00:16:39.990000Z | 100.0 Hz, 100000 samples\n"
            "XX.STA.00.BX2 | 2020-01-01T00:00:00.000000Z - 2020-01-01T00:16:39.990000Z | 100.0 Hz, 100000 samples\n"
            "XX.STA.00.BX3 | 2020-01-01T00:00:00.000000Z - 2020-01-01T00:16:39.990000Z | 100.0 Hz, 100000 samples"
        )
        str_cleaned = (
            "4 Trace(s) in Stream:\n"
            "XX.STA.00.BDH       | 2020-01-01T00:00:00.000000Z - 2020-01-01T00:16:39.990000Z | 100.0 Hz, 100000 samples\n"
            "XX.STA.00-H.BX1     | 2020-01-01T00:00:00.000000Z - 2020-01-01T00:16:39.990000Z | 100.0 Hz, 100000 samples\n"
            "XX.STA.00-H-1.BX2   | 2020-01-01T00:00:00.000000Z - 2020-01-01T00:16:39.990000Z | 100.0 Hz, 100000 samples\n"
            "XX.STA.00-H-1-2.BX3 | 2020-01-01T00:00:00.000000Z - 2020-01-01T00:16:39.990000Z | 100.0 Hz, 100000 samples"
        )
        self.assertEqual(self.stream.__str__(), str_standard)
        obj = CleanedStream(self.stream)
        assert isinstance(obj, CleanedStream)
        assert isinstance(obj, Stream)
        self.assertEqual(obj.__str__(), str_cleaned)
        # Test lossless back-and-forth
        obj = Stream(obj)
        assert not isinstance(obj, CleanedStream)
        assert isinstance(obj, Stream)
        self.assertEqual(obj.__str__(), str_standard)
        obj = CleanedStream(obj)
        assert isinstance(obj, CleanedStream)
        assert isinstance(obj, Stream)
        self.assertEqual(obj.__str__(), str_cleaned)

    def test_tag(self):
        """Test wrapping of tag() method"""
        obj = CleanedStream(self.stream)
        obj = obj.tag('ROT', cleaned_ids=('*BX3', ))
        self.assertEqual(
            obj.__str__(),
            "4 Trace(s) in Stream:\n"
            "XX.STA.00.BDH           | 2020-01-01T00:00:00.000000Z - 2020-01-01T00:16:39.990000Z | 100.0 Hz, 100000 samples\n"
            "XX.STA.00-H.BX1         | 2020-01-01T00:00:00.000000Z - 2020-01-01T00:16:39.990000Z | 100.0 Hz, 100000 samples\n"
            "XX.STA.00-H-1.BX2       | 2020-01-01T00:00:00.000000Z - 2020-01-01T00:16:39.990000Z | 100.0 Hz, 100000 samples\n"
            "XX.STA.00-H-1-2-ROT.BX3 | 2020-01-01T00:00:00.000000Z - 2020-01-01T00:16:39.990000Z | 100.0 Hz, 100000 samples"
        )

    # def test_plot(self):
    #     """Just makes sure that nothing blows up"""
    #     obj = CleanedStream(self.stream)
    #     obj.plot()

def suite():
    return unittest.makeSuite(TestMethods, "test")


if __name__ == "__main__":
    unittest.main(defaultTest="suite")
