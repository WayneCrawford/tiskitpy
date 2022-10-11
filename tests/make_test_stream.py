#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Test stream for the test_* routines
"""
from dataclasses import dataclass
import random
from copy import copy

import numpy as np
from obspy.core import UTCDateTime
from obspy.core.stream import Stream, Trace

@dataclass
class SineParms:
    frequency: float
    amplitude: float
    phase: float
    
    def make_sine(self, t):
        """Make a sine wave
        
        Args:
            t: time array
        """
        return self.amplitude * np.sin(2*np.pi*self.frequency*t + np.radians(self.phase))

    def make_sine(self, t):
        """Make a sine wave
        
        Args:
            t: time array
        """
        return self.amplitude * np.sin(2*np.pi*self.frequency*t + np.radians(self.phase))

# default_sineparms = {
#     "BX1": [SineParms(13., 10., 0), SineParms(10., 0., 0.), SineParms(42., 0., 0.)],
#     "BX2": [SineParms(13., 2., 45.), SineParms(10., 1., 90.), SineParms(42., 3., 135.)],
#     "BX3": [SineParms(13., 10., 0), SineParms(10., 10., 0), SineParms(42., 10., 0)]}
sineparms = {
    "tilt": [SineParms(13., .33, 0), SineParms(10., .33, 0.), SineParms(42., .33, 0.)],
    "compliance": [SineParms(7, .33, 45.), SineParms(23, .33, 90.), SineParms(31, .34, 135.)],
    "other":      [SineParms(0, 1., 0)]
}

chans = {"BX1": {'tilt': 1, 'compliance': 1e-9, 'other': .01},
         "BX2": {'tilt': .9, 'compliance': 1e-9, 'other': .01},
         "BX3": {'tilt': .1, 'compliance': .1,   'other': .03},
         "BDH": {'tilt': 1e-4, 'compliance': 1,   'other': .01}}
                 

def make_test_stream(sineparms=sineparms,
                     net="XX", sta="STA", data_len=1000.0, sampling_rate=100,
                     loc="00", starttime=UTCDateTime(2020, 1, 1)):
    """
    Sum of sinusoids with different amplitudes and phases
    
    Returns:
        (tuple):
            stream (:class:`obspy.stream.Stream`): stream with channels
              as given in chans
            sp (dict): for each stream channnel, an array of the constituant
                SineParms parameters
    """
    
    stream = Stream()
    header = dict(sampling_rate=sampling_rate, network=net, station=sta,
                  location=loc, starttime=starttime)
    sp = {}
    t = np.arange(0, data_len, 1.0 / sampling_rate)
    for chan in chans.keys():
        ch_sp = []
        header["channel"] = chan
        data = np.zeros(t.shape)
        for key, ampl in chans[chan].items():
            for s in sineparms[key]:
                ch_s = copy(s)
                ch_s.amplitude *= ampl
                if ch_s.frequency == 0:
                    ch_s.frequency = random.random()*sampling_rate/2.
                data += ch_s.make_sine(t)
                ch_sp.append(ch_s)
        trace = Trace(data, header)
        stream += trace
        sp[trace.id] = ch_sp
    return stream, sp

