from tiskitpy.rptransient import PeriodicTransient
from obspy.core.stream import read


transient_starttime = '2009-01-23T09:00:00'
transient_clips = (-1000, 1000)
transient_period = 3600
dp = 1

stream = read('data/LSVSB.Z.sample.mseed', 'MSEED')
trace_train = stream.select(channel="*MHZ")[0].copy()

# PERIODIC TRANSIENT
pt = PeriodicTransient("hourly_glitch",
                       transient_period,
                       dp=dp, clips=transient_clips,
                       transient_starttime=transient_starttime)

# Calculate timing
pt.calc_timing(trace_train)

# Calculate the transient
pt.calc_transient(trace_train, plots=True)

# Remove the transient
stream_before = stream.copy()
trace = stream.select(channel="*MHZ")[0]
trace_after = pt.remove_transient(trace, plots=True)

# Replace the Z component trace with the one with the removed transient
for i, t in enumerate(stream):
    if t.stats.channel == 'MHZ':
        print('replacing channel')
        stream[i] = trace_after
        break

# Compare filtered stream before and after
trace_before = stream_before.select(channel='MHZ').copy()
trace_after = stream.select(channel='MHZ').copy()
trace_before[0].stats.location='99'
stream_compare = trace_before + trace_after
stream_compare.filter("lowpass", freq=pt.freq_LP)
stream_compare.plot()

