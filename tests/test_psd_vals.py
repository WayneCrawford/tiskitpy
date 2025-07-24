#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Functions to test the lcheapo functions
"""
# from os import system
import sys
import unittest
import filecmp
import inspect
import difflib
from pathlib import Path
import pytest
import pickle
import datetime
from copy import deepcopy

from obspy.core.inventory.response import Response
from obspy import read_inventory, UTCDateTime
from obspy.core.stream import read as stream_read, Stream
from matplotlib import pyplot as plt
import numpy as np

from tiskitpy.synthetic import PSDVals
from tiskitpy import SpectralDensity


def _sloped_psd_vals(val_1Hz, slope, logf_low, logf_high, logf_step, units="m/s^2"):
    # Make list of [f, level] for a linear slope
    x = [[f, val_1Hz * np.power(f, slope)]
         for f in np.power(10, np.arange(logf_low, logf_high, logf_step))]
    # Create and return PSDVals object
    return PSDVals((x, False), units)


@pytest.mark.parametrize(
    "psdv",
    [PSDVals(PSDVals.sloped_freqs_and_values(-100, 20, -4, 0.1, 1), 'm/s^2'),
     PSDVals(PSDVals.sloped_freqs_and_values(-120,   0, -4, 0.1, 1), 'm/s^2'),
     PSDVals(PSDVals.sloped_freqs_and_values(-140, -20, -4, 0.1, 1), 'm/s^2')])
def test_accel_as_vel(psdv):
    """
    Verify that converting from accel to vel has expected behavior
    """
    plotit = False

    stats = {'starttime': UTCDateTime(2010, 1, 1),
             'endtime':   UTCDateTime(2010, 1, 10),
             'sampling_rate': 1.}
    # Blue acceleration spectrum
    # psdv = _sloped_psd_vals(2*np.pi*10.**-6, 1, -4, 0.1, 1)
    
    if plotit is True:
        # Plot the normal response, and the .accel_as_vel response
        fig, ax = plt.subplots()
        ax.semilogx(psdv.freqs, psdv.values, c='b', label='as accel')
        ax.semilogx(psdv.accel_as_vel.freqs, psdv.accel_as_vel.values,
                c='r', label='accel_as_vel')
        ax.axvline(1/(2*np.pi), ls='--')
        ax.set_xlabel('Frequency')
        ax.set_ylabel('dB ref UNITS')
        plt.legend()
        plt.show()

    # transform to velocity and use response w.r.t. vel
    vel_gain = 1500.
    stats['response'] = Response.from_paz([0, 0], [0.01, 0.01, 20, 20], vel_gain,
                                          input_units='m/s', output_units='count')
    if plotit is True:
        stats['response'].plot
    trace_vel, phases = psdv.accel_as_vel.as_trace(stats, channel='LHZ', location='01')
    # Keep as accel and use response w.r.t. accel
    stats['response'] = Response.from_paz([0], [0.01, 0.01, 20, 20], vel_gain,
                                          input_units='m/s**2', output_units='count')
    if plotit is True:
        stats['response'].plot
    # Need to add pi/2 to vel phases for accel phases
    trace_accel, _ = psdv.as_trace(stats, channel='LHZ', location='02', phases=phases+np.pi/2)

    sd = SpectralDensity.from_stream(Stream([trace_vel, trace_accel]),
                                     windowtype='prol4pi')
    # Compare the traces made with each method
    if plotit is True:
        fig, ax = plt.subplots()
        ax.plot(trace_vel.times()/60, trace_vel.data, 'b-', label='vel_based (loc 01)')
        ax.plot(trace_accel.times()/60,trace_accel.data, 'r:', label='accel_based (loc 02)')
        ax.set_xlabel('Minutes')
        ax.set_xlim((0, 5))
        ax.set_ylabel('Velocity Amplitude')
        plt.legend()
        plt.show()
        ax = sd.plot(show=False)
        _fix_and_plot(ax[0, 0][0], psdv.freqs, psdv.values, '+')
        _fix_and_plot(ax[0, 1][0], psdv.freqs, psdv.values, '+')
        plt.show()
        
    # Test that spectra are similar
    threshold = 0.05  # acceptable ratio difference between accel-based and vel-based spectra
    # Spectra of time series created from acceleration and from velocity
    spect_diff = np.abs(sd.autospect('*.01.LHZ'))-np.abs(sd.autospect('*.02.LHZ'))
    spect_sum = np.abs(sd.autospect('*.01.LHZ'))+np.abs(sd.autospect('*.02.LHZ'))
    spect_ratio = np.abs(spect_diff / spect_sum)
    assert np.mean(np.abs(spect_ratio)) <  threshold
    # Spectra of time series versus input model
    model_interp = np.interp(np.log10(sd.freqs), np.log10(psdv.freqs), psdv.values)
    model_interp = np.power(10., model_interp/10.)
    for seed_id in ('*.01.LHZ', '*.02.LHZ'):
        spect_diff = np.abs(sd.autospect(seed_id)) - model_interp
        spect_ratio = np.abs(spect_diff / model_interp)
        assert np.mean(np.abs(spect_ratio)) <  threshold


def _fix_and_plot(ax, xs, ys, sym):
    xlim, ylim = ax.get_xlim(), ax.get_ylim()
    ax.plot(xs, ys, sym)
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)


class TestMethods(unittest.TestCase):
    """
    Test suite
    """
    def setUp(self):
        self.path = Path(inspect.getfile(
            inspect.currentframe())).resolve().parent
        self.test_path = self.path / "data" / "decimate"

    def test_add_sub(self):
        obj = PSDVals(([[0.001, -100], [1, -10]], True), 'm/s^2')
        self.assertListEqual(list(obj.freqs), [0.001, 1])
        self.assertListEqual(list(obj.values), [-100, -10])
        # __add()__
        obj += 20
        self.assertListEqual(list(obj.values), [-80, 10])
        # __sub()__
        obj -= 40
        self.assertListEqual(list(obj.values), [-120, -30])
        # accel_as vel (property)
        self.assertListEqual(list(obj.accel_as_vel.values),
                             list(np.array([-120, -30]) -
                                  20*np.log10(2*np.pi*np.array([0.001, 1]))))

    def test_resample(self):
        obj = PSDVals(([[0.001, -120], [1, -30]], True), 'm/s^2')
        print(obj)
        resamp_freqs = np.array([0.001, 0.01, 0.1, 1])
        resamp_expected_values = [-120, -90, -60, -30]

        # resample_values()
        new_freqs = np.array([0.001, 0.01, 0.1, 1])
        self.assertListEqual(list(obj.resample_values(resamp_freqs)),
                             resamp_expected_values)
        # resample()
        obj.resample(resamp_freqs)
        self.assertListEqual(list(obj.values), resamp_expected_values)

    def test_as_trace_spectra_no_response(self):
        """
        Verify that created time series has expected spectral shape
        
        Spectra with no instrument response
        There are no assertions, just a possible visual vailidation
        """
        plotit = False # Visual validation of what's going on
        
        # "Red" spectra
        traces = {}
        psds = {'normal': _sloped_psd_vals(10., -1.5, -3, 0.1, 0.25),
                'wide':   _sloped_psd_vals(10., -1.5, -4, 1.0, 0.25),
                'spaced': _sloped_psd_vals(10., -1.5, -3, 0.1, 1)}
        for k, obj in psds.items():
            traces[k], _ = obj.as_trace(dict(starttime=UTCDateTime(2010,1,1),
                                             endtime=UTCDateTime(2010,1,5),
                                             sampling_rate=1.),
                                     plotit=plotit)
            sd = SpectralDensity.from_stream(Stream([traces[k]]), windowtype='prol4pi')
        if plotit is True:
            _plot_all(traces, psds, 'red')
        
        # "Blue" spectra
        traces = {}
        psds = {'normal': _sloped_psd_vals(10., 1.5, -3, 0.1, 0.25),
                'wide':   _sloped_psd_vals(10., 1.5, -4, 1.0, 0.25),
                'spaced': _sloped_psd_vals(10., 1.5, -3, 0.1, 1)}
        for k, obj in psds.items():
            traces[k], _ = obj.as_trace(dict(starttime=UTCDateTime(2010,1,1),
                                             endtime=UTCDateTime(2010,1,5),
                                             sampling_rate=1.),
                                     plotit=plotit)
            sd = SpectralDensity.from_stream(Stream([traces[k]]), windowtype='prol4pi')
        if plotit is True:
            _plot_all(traces, psds, 'blue')

        # "White" spectra
        traces = {}
        psds = {'normal': _sloped_psd_vals(10., 0., -3, 0.1, 0.25),
                'wide':   _sloped_psd_vals(10., 0., -4, 1.0, 0.25),
                'spaced': _sloped_psd_vals(10., 0., -3, 0.1, 1)}
        for k, obj in psds.items():
            traces[k], _ = obj.as_trace(dict(starttime=UTCDateTime(2010,1,1),
                                             endtime=UTCDateTime(2010,1,5),
                                             sampling_rate=1.),
                                     plotit=plotit)
            sd = SpectralDensity.from_stream(Stream([traces[k]]), windowtype='prol4pi')
        if plotit is True:
            _plot_all(traces, psds, 'white')
       

def suite():
    return unittest.makeSuite(TestMethods, 'test')


def _plot_all(traces, psds, shape_text):
    f, ax = plt.subplots()
    for k, tr in traces.items():
        sd_4pi = SpectralDensity.from_stream(Stream([tr]), windowtype='prol4pi')
        sd_1pi = SpectralDensity.from_stream(Stream([tr]), windowtype='prol1pi')
        lines = ax.semilogx(sd_1pi.freqs, 10*np.log10(sd_1pi.autospect(sd_1pi.ids[0])), ':', label=k+'_out_1pi')
        c = lines[0].get_color()
        ax.semilogx(sd_4pi.freqs, 10*np.log10(sd_4pi.autospect(sd_4pi.ids[0])), c=c, label=k+'_out_4pi')
        ax.semilogx(psds[k].freqs, psds[k].values, c=c, ls='None', marker='+', label=k+'_input')
    plt.suptitle(shape_text)
    ax.legend()
    plt.show()


if __name__ == '__main__':
    unittest.main(defaultTest='suite')
