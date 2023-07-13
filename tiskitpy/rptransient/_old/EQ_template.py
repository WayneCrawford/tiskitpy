#!env python3
"""Create a template of time ranges to skip after earthquakes

REPLACED BY EQRemover
"""

#################################################
# import obspy
from obspy.clients.fdsn import Client
from obspy.core import UTCDateTime
from obspy.core.stream import Stream, Trace
import numpy as np
import matplotlib.pyplot as plt


class EQTemplate(Trace):
    """
    A template for zeroing out data after big earthquakes
    """
    def __init__(self, inp, start_time, end_time, min_magnitude=5.85,
                 days_per_magnitude=1.5,
                 plot=False, verbose=True, eqfile=None):
        """
        Make a template to cut out times after large EQs

        Arguments:
            inp: Input data (opspy.core.Trace)
            start_time (obspy.core.UTCDateTime): earliest data that will be
                presented
            end_time (obspy.core.UTCDateTime): latest data that will be
                presented
            min_magnitude (float): EQ Magnitude above which to cut out times
            days_per_magnitude (float): days to cut per magnitude above
                min_magnitude
            verbose: print out list of EQs and the time that will be cut
            eqfile: name of  csv file containing global earthquakes
                (column 1 = time (ISO8601), column4=magnitude)
                If this is not provided, the program will download the data
                based on starttime and endtime.  It will load early enough
                before to cover an M9 EQ
        :returns: template (obspy.core.Trace)
        
        Here's an obspy example of downloading the eqfile:
        from obspy import UTCDateTime
        """
        cat = Client("USGS").get_events(
            starttime=start_time-86400*(9-min_magnitude)*days_per_magnitude,
            endtime=end_time,
            minmagnitude=min_magnitude, 
            orderby='time-asc')
        f = open(eqfile, 'r')
        for line in f:
            w = line.split(',')
            if w[0] == 'time':
                continue
            starttime = UTCDateTime(w[0])
            mag = float(w[4])
            if starttime >= inp.stats.endtime:
                continue
            cutdays = (mag-min_magnitude) * days_per_magnitude
            # This is done before the magnitude limit check simply to allow
            #    printout of EQs within the acceptable time frame
            if cutdays < 0:
                cutdays = 0.01
            endtime = starttime + cutdays * 86400
            if endtime <= self.stats.starttime:
                continue
            print('M{:.1f} EQ of {} ...'.format(mag, w[0]), end='')
            if mag <= min_magnitude:
                print('magnitude too small to worry about')
                continue
            print('blanking {:.1f} days'.format(cutdays))
            startaddr = int((starttime - self.stats.starttime) * sps)
            if startaddr < 0:
                startaddr = 0
            if endtime > self.stats.endtime:
                endtime = self.stats.endtime
            endaddr = int((endtime - inp.stats.starttime) * sps)
            self.data[startaddr:endaddr] = 0
        if plot:
            self.plot_mask(inp)
            
    def apply(self, inp, plot=False):
        """
        Zero out data in trace when it is too close to an eqfile
        
        Replaced by EQRemover.zero()

        Returns an error if the EQ catalog does not bound the Trace time
        Returns a warning if the EQ catalog does not start far enough before
        the trace to cover an M9 EQ.
        
        Arguments:
            inp (obspy.core.Trace): seismological data
            plot: plot traces with EQs cut out
        """
        if not isinstance(inp, Trace):
            raise(ValueError, 'inp is not an obspy Trace object')
        super().__init__(np.ones(inp.data.shape), inp.stats)
        self.stats.channel = 'TTT'
        sps = inp.stats.sampling_rate

         

    def __mul__(self, val):
        """
        Return data multiplied by template
        
        Replaced by EQRemover.zero()
        
        :type val: Stream or Trace
        """
        assert isinstance(val, (Trace, Stream)), f'{type(val)=}'
        newval = val.copy()
        if isinstance(val, Trace):
            assert self.data.size == val.data.size
            newval.data = self.data * val.data
        elif isinstance(val, Stream):
            for oldtrace, newtrace in zip(val, newval):
                assert self.data.size == oldtrace.data.size
                newtrace.data =  self.data * oldtrace.data
        return val

    def plot_mask(self, inp):
        """
        Plot template times inp
        """
        modf = inp.copy()
        modf.data = self.data * inp.data
        modf.stats.channel = 'XX' + inp.stats.channel[2]
        fig = plt.figure(num='EQ_template')
        Stream([modf, inp, self]).plot(color='blue', fig=fig, equal_scale=False)
        fig.show()


#########################################################################
# if __name__ == "__main__":
# 	sys.exit(main())
