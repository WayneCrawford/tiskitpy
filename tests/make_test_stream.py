#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Test stream for the test_* routines
"""
from dataclasses import dataclass

import numpy as np
from obspy.core import UTCDateTime
from obspy.core.stream import Stream, Trace

@dataclass
class SineParms:
    frequency: float
    amplitude: float
    phase: float

default_sineparms = {"BX1": [SineParms(13., 10., 0),
                             SineParms(10., 0., 0.),
                             SineParms(42., 0., 0.)],
                     "BX2": [SineParms(13., 2., 45.),
                             SineParms(10., 1., 90.),
                             SineParms(42., 3., 135.)],
                     "BX3": [SineParms(13., 10., 0),
                             SineParms(10., 10., 0),
                             SineParms(42., 10., 0)]}


def make_test_stream(sineparms=default_sineparms,
                     net="XX", sta="STA", data_len=1000.0, sampling_rate=100,
                     loc="00", starttime=UTCDateTime(2020, 1, 1)):
    """
    Sum of sinusoids with different amplitudes and phases
    """
    stream = Stream()
    header = dict(sampling_rate=sampling_rate, network=net, station=sta,
                  location=loc, starttime=starttime)
    sp = {}
    for chan in sineparms.keys():
        header["channel"] = chan
        t = np.arange(0, data_len, 1.0 / sampling_rate)
        data = np.zeros(t.shape)
        for s in sineparms[chan]:
            data += s.amplitude * np.sin(2*np.pi*s.frequency*t
                                         + np.radians(s.phase))
        trace = Trace(data, header)
        stream += trace
        sp[trace.id] = sineparms[chan]
    return stream, sp

