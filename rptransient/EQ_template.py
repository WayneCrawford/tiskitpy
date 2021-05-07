#!env python3
"""Create a template of time ranges to skip after earthquakes"""

#################################################
# import obspy
from obspy.core import UTCDateTime
from obspy.core.stream import Stream, Trace
import numpy as np
import matplotlib.pyplot as plt

def_mag_limit = 5.85
def_days_per_magnitude = 1.5


class EQTemplate(Trace):
    """
    An obspy Trace with zeros after big EQs, ones elsewhere
    """
    def __init__(self, inp, eqfile, mag_limit=def_mag_limit,
                 days_per_magnitude=def_days_per_magnitude,
                 plot=False, verbose=True):
        """
        Make a template to cut out times after large EQs

        :param inp: seismological channel (obspy.core.Trace)
        :param eqfile: name of  csv file containing global earthquakes
            (column 1 = time (ISO8601), column4=magnitude)
            Can be downloaded from USGS website or using web services 
            obspy example:
                obspy.clients.fdsn.client.Client("IRIS").get_events())
        :param mag_limit: EQ Magnitude above which to cut out times
        :param days_per_magnitude: days to cut per magnitude above mag_limit
        :param verbose: print out list of EQs and the time that will be cut
        :param plot: plot traces with EQ cut out
        :returns: template (obspy.core.Trace)
        
        Here's an obspy example of downloading the eqfile:
        from obspy.clients.fdsn import Client
        from obspy import UTCDateTime
        Client("USGS").get_events(starttime=UTCDateTime('2016-02-27'),
                                  endtime=UTCDateTime('2017-02-28'),
                                  minmagnitude=5.5, format='csv', orderby='time-asc',
                                  filename='eqs_20160227-20170228.csv')
        """
        assert isinstance(inp, Trace)
        super().__init__(np.ones(inp.data.shape), inp.stats)
        self.stats.channel = 'TTT'
        sps = inp.stats.sampling_rate

        f = open(eqfile, 'r')
        for line in f:
            w = line.split(',')
            if w[0] == 'time':
                continue
            starttime = UTCDateTime(w[0])
            mag = float(w[4])
            if starttime >= inp.stats.endtime:
                continue
            cutdays = (mag-mag_limit) * days_per_magnitude
            # This is done before the magnitude limit check simply to allow
            #    printout of EQs within the acceptable time frame
            if cutdays < 0:
                cutdays = 0.01
            endtime = starttime + cutdays * 86400
            if endtime <= self.stats.starttime:
                continue
            print('M{:.1f} EQ of {} ...'.format(mag, w[0]), end='')
            if mag <= mag_limit:
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
            
    def __mul__(self, val):
        """
        Return data multiplied by template
        
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
