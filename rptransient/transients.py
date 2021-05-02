#!env python3
""" Model and remove glitches from BBOBS data"""

#################################################
# import obspy
from obspy.core import UTCDateTime
from obspy.core.stream import Stream
import diracComb_v2 as dc
import matplotlib.pyplot as plt
import numpy as np
# from scipy import signal
# import math as M
# import sys
# import glob

from .periodic_transient import PeriodicTransient
from .utilis import prep_filter

def_mag_limit = 5.85
def_days_per_magnitude = 1.5

class Transients():
    """
    Periodic Transients class
    """
    def __init__(trans_parms=[]]):
        """
        :param trans_parms: list of transient parameters
        """
        self.trans_parms = trans_parms

    def calc_times(trace, eqfile=None, plots=False, mag_limit=def_mag_limit,
                   days_per_magnitude=def_days_per_magnitude):
        """
        Calculate transient times

        :param trace: data
        :param eqfile: earthquake catalog file (to blank noisy times)
                       should probably be changed to obspy catalog
        """
        starttime = trace.stats.starttime
        endtime = trace.stats.endtime

        # sps = trace.stats.sampling_rate
        station = trace.stats.station
        print('Station='+station)

        #######################################################################
        # WEILANDT METHOD TO REMOVE GLITCHES
        ######################################################################
        filt_trace = prep_filter(trace)

        # Make a template to cut out EQ times, clip extreme values...)
        template = _make_EQ_template(filt_trace, eqfile, mag_limit,
                                     days_per_magnitude, plotit=plots)

        # Loop through different glitch periods (1 hour, 3 hours, ...) and
        # remove the glitches
        iGlitchType = 0
        glitches, dirac_combs, nGlitcheses = [], [], []
        for gp in self.trans_parms:
            iGlitchType += 1
            print('='*75)
            print(' '*15 + 'Fitting {:g}-second glitch'.format(gp.period))
            print('='*75)

            ref_time, shifts = self._get_transient_shifts(gp, filt_trace)

            if shifts:
                starttimes = shifts.copy().insert(0, starttime)
                endtimes = shifts.copy().append(endtime)
                reftimes = shifts.copy().insert(0, ref_time)
                print('Should have some way to cut out last period (because will '
                      'be full of leveling)')
                for s, e, r in zip(starttimes, endtimes, reftimes):
                    tr_slice = trace.cut(starttime=s, endtime=e)
                    te_slice = template.cut(starttime=s, endtime=e)
                    gp.verify(tr_slice, te_slice, r)
            else:
                gp.verify(trace, template, ref_time)

        return self.trans_parms
        
    def calc_trans(stream, times=[], catalog=None):
        """
        Return transient model

        :param stream: data (trace?)
        :param times: list of transient times
        :param catalog: earthquake catalog (to blank noisy times) 
        """
        """
        Calculate transient times

        :param trace: data
        :param eqfile: earthquake catalog file (to blank noisy times)
                       should probably be changed to obspy catalog
        """
        starttime = trace.stats.starttime
        endtime = trace.stats.endtime

        # sps = trace.stats.sampling_rate
        station = trace.stats.station
        print('Station='+station)

        #######################################################################
        # WEILANDT METHOD TO REMOVE GLITCHES
        ######################################################################
        filt_trace = prep_filter(trace)

        # Make a template to cut out EQ times, clip extreme values...)
        template = _make_EQ_template(filt_trace, eqfile, mag_limit,
                                     days_per_magnitude, plotit=plots)

        # Loop through different glitch periods (1 hour, 3 hours, ...) and
        # remove the glitches
        for gp in self.trans_parms:
            print('='*75)
            print(' '*15 + 'Fitting {gp}')
            print('='*75)

            ref_time, shifts = gp._get_transient_shifts(filt_trace)

            if shifts:
                starttimes = shifts.copy().insert(0, starttime)
                endtimes = shifts.copy().append(endtime)
                reftimes = shifts.copy().insert(0, ref_time)
                print('Should cut last period (will be full of leveling)')
                for s, e, r in zip(starttimes, endtimes, reftimes):
                    tr_slice = trace.cut(starttime=s, endtime=e)
                    te_slice = template.cut(starttime=s, endtime=e)
                    gp.verify(tr_slice, te_slice, r)
            else:
                gp.verify(trace, template, ref_time)

        return self.trans_parms
        
    def clean_trans(stream, trans_model, times):
        """
        Remove transient from data

        :param stream: data (trace?)
        :param trans_model: transient_model
        :param times: list of transient times
        :param catalog: earthquake catalog (to blank noisy times) 
        """


def removeGlitch(stream, trans_parms, eqfile, interactive=False, plots=False,
                 glitch=None, diracComb=None, match=True,
                 mag_limit=def_mag_limit, days_per_magnitude=def_days_per_magnitude):
    """
    Calculate and remove periodic glitches from one seismo channel

    :param trace: seismological channel
    :param trans_parms = list of TransParms objects
    :param eqfile = csv file containing global earthquakes
                (column 1 = time (ISO8601), column4=magnitude)
    :param interactive (boolean): whether to let the user view and modify
                                  glitch interval and clip levels
    :param plots (boolean):  Whether to plot stuff
    #### THESE NEXT TWO ARE NOT YET IMPLEMENTED, MEANT FOR REMOVE
    #### GLITCHES BASED ON PRIOR INFORMATION
    :param glitches: glitches to use (will calculate the best multiplication
                     factor to fit the given trace) [None]
    :type glitches: list of numpy.ndarrays
    :param diracComb (list of obspy traces): traces to convolve with glitches
                to generate synthetic glitch time series [None]
    :param match: Match and cancel each pulse separately
    :type match: bool

    Outputs:
        new_stream : (obspy.trace) glitch-corrected trace
        glitches   : (list of numpy.ndarray) list of average Glitches (one per
                     element in trans_parms)
        diracCombs : (list of obspy traces): forms to convolve with
                      glitch to get synthetic glitch time series
        nGlitcheses: (list of ints): number of data glitches used to calculate
                                    average glitch(es)
    """
    starttime = stream[0].stats.starttime
    endtime = stream[0].stats.endtime

    trace = stream.select(component='Z')[0]
    # sps = trace.stats.sampling_rate
    station = trace.stats.station
    print('Station='+station)

    #######################################################################
    # WEILANDT METHOD TO REMOVE GLITCHES
    ######################################################################
    # Bandpass to 2-30 mHz: demean; LP 30 mHz, order 6 (dt=1.6s);
    # HP 2 mHz, order 2 (dt=1.6s)
    trace.detrend('demean')
    trace.detrend('linear')
    trace = trace.filter('lowpass', freq=0.03, corners=6)
    trace = trace.filter('highpass', freq=0.002, corners=2)

    # Make a template to cut out EQ times, clip extreme values...)
    template = _make_EQ_template(trace, eqfile, mag_limit, days_per_magnitude, plotit=plots)

    # Loop through different glitch periods (1 hour, 3 hours, ...) and
    # remove the glitches
    iGlitchType = 0
    glitches, dirac_combs, nGlitcheses = [], [], []
    for gp in trans_parms:
        iGlitchType += 1
        print('='*75)
        print(' '*15 + 'Fitting {:g}-second glitch'.format(gp.period))
        print('='*75)

        shifts = []
        if isinstance(gp.starttime, list):
            # Search for glitch start corresponding to data start and
            # make a list of glitch start shifts within data
            ref_time = gp.starttime
            if len(gp.starttime) > 1:
                for test_time in gp.starttime[1:]:
                    if test_time < starttime:
                        ref_time = test_time
                    elif test_time < endtime:
                        # there is a break in glitch times within this data
                        shifts.append(test_time)
        elif gp.starttime:
            ref_time = gp.starttime
        else:
            ref_time = trace.stats.starttime

        if shifts:
            starttimes = shifts.copy().insert(0, starttime)
            endtimes = shifts.copy().append(endtime)
            reftimes = shifts.copy().insert(0, ref_time)
            print('Should have some way to cut out last period (because will '
                  'be full of leveling)')
            for s, e, r in zip(starttimes, endtimes, reftimes):
                tr_slice = trace.cut(starttime=s, endtime=e)
                te_slice = template.cut(starttime=s, endtime=e)
                out, gls, dcs, ngs = _calcglitch(tr_slice, te_slice, r, gp,
                                                 iGlitchType, match,
                                                 interactive, plots)
        else:
            out, gls, dcs, nGs = _calcglitch(trace, template, ref_time, gp,
                                             iGlitchType, match,
                                             interactive, plots)
        nGlitcheses.append(nGs)
        glitches.append(gls)
        dirac_combs.append(dcs)

    return (out, glitches, dirac_combs, nGlitcheses)


def plot_EQ_template(trace, eqfile, mag_limit=def_mag_limit, days_per_magnitude=def_days_per_magnitude,
                     verbose=True):
    template = _make_EQ_template(trace, eqfile, mag_limit, days_per_magnitude, verbose)
    temp = template.copy()
    temp.data = temp.data * trace.data
    Stream([trace, temp]).plot(color='blue')

def _make_EQ_template(trace, eqfile, mag_limit=def_mag_limit, days_per_magnitude=def_days_per_magnitude,
                      verbose=True):
    """
    Make a template to cut out times after large EQs

    :param trace: seismological channel (obspy.trace)
    :param eqfile: name of  csv file containing global earthquakes
                   (column 1 = time (ISO8601), column4=magnitude)
    :param mag_limit: EQ Magnitude above which to cut out times
    :param days_per_magnitude: time (in days) to cut for each magnitude above mag_limit
    :param verbose: print out list of EQs and the time that will be cut [True]
    :param plotit: plot traces with EQ cut out [True]

    :returns: template (obspy.trace) of ones for data to keep, zeros for
              data to remove
    """
    template = trace.copy()
    template.data = np.ones(trace.data.shape)
    template.stats.channel = 'TTT'
    sps = trace.stats.sampling_rate

    f = open(eqfile, 'r')
    for line in f:
        w = line.split(',')
        if w[0] == 'time':
            continue
        starttime = UTCDateTime(w[0])
        mag = float(w[4])
        if starttime >= trace.stats.endtime:
            continue
        cutdays = (mag-mag_limit) * days_per_magnitude
        # This is done before the magnitude limit check simply to allow
        #    printout of EQs within the acceptable time frame
        if cutdays < 0:
            cutdays = 0.01
        endtime = starttime + cutdays * 86400
        if endtime <= trace.stats.starttime:
            continue
        print('M{:.1f} EQ of {} ...'.format(mag, w[0]), end='')
        if mag <= mag_limit:
            print('magnitude too small to worry about')
            continue
        print('blanking {:.1f} days'.format(cutdays))
        startaddr = int((starttime - trace.stats.starttime) * sps)
        if startaddr < 0:
            startaddr = 0
        if endtime > trace.stats.endtime:
            endtime = trace.stats.endtime
        endaddr = int((endtime - trace.stats.starttime) * sps)
        template.data[startaddr:endaddr] = 0
    return template


#########################################################################
# if __name__ == "__main__":
# 	sys.exit(main())
