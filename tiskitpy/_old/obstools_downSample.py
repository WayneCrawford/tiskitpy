#!env python3
""" Resample data to a easy-to-handle rate"""

###############################################################################
# Note: resample takes too long, even on 1.25 sps data
#       Could (should?) write:
#         -  resample extension that uses FIR filter (min phase or zero phase,
#            use scipy.signal.firwin to calculate FIR)
#         - decimate extension that uses FIR filter (min or zero phase,
#           scipy.signal.firwin)
#       Both should return the FIR filter, which could be stuffed into the
#       instruments' response information
#       scipy doesn't seem to have a minimum phase FIR method: could use Scherbaum
#       method to  convert a zero phase to minimum phase
###############################################################################
import obspy
from obspy.clients.filesystem import sds
import time
import glob

###############################################################################
def decimate(stream,decimates,verbose=False):
    """
    Decimate obspy stream in some kind of intelligent way
    
    Usage: stream_decim = decimate(stream,decimates)
    
    Inputs:
        stream (obspy.stream)
        decimates: list of decimation factors (use 2,3 or 5)
        
    Outputs:
        stream_decim (obspy.stream) decimated data
    """
    sr=stream[0].stats.sampling_rate
    decimate_total=1
    for d in decimates: 
        decimate_total=decimate_total*d

    print('Decimating data from {:g} to {:g} Hz ({:d}x)... '.format(sr,
                                sr/decimate_total, decimate_total),end='',flush=True)
    tic=time.time()
    # This method from Martha Deen, except I removed 'no_filter', which is  faster
    # but aliased!!!!!  Was she was counting on final (very slow) resample() to take
    # care of filtering???????!?!?!?
    for d in decimates:
        stream.decimate(d)
        if verbose:
            print(stream)
            for trace in stream:
                print(trace)
                print(trace.stats.processing)
            
    print('Took {:.1f} seconds'.format(time.time()-tic))
    print('New data has {} samples'.format([tr.data.size for tr in stream]))
    stream.verify()
    return stream
    
###############################################################################
def readMSEED(filenames,starttime,endtime,verbose=False):
    """
    Read in miniSEED files, using either obspy.read or, if the file is 
    more than 2.1 Gb (libmseed limitation), by reading in chunks first
    """
    import io
    import os.path
    
    firstTime=True
    for filename in glob.glob(filenames) :
        filesize=os.path.getsize(filename)
        if filesize < 2^31 :
            if verbose:
                print("Filesize {:.2f}Go < 2^31 bytes, reading using obspy.read".format(\
                        filesize/float(1024*1024*1024)))
            stream= obspy.read(filename,format='MSEED', starttime=starttime,endtime=endtime)
        else:
            if verbose:
                print("Filesize {:.2f}Go >= 2^31 bytes, breaking into readable chunks for obspy.read".format(\
                        filesize/float(1024*1024*1024)))
            reclen = 4096
            chunksize = 20000 * reclen # Around 80 MB
            madeStream=False

            with io.open(filename, "rb") as fh:
                while True:
                    with io.BytesIO() as buf:
                        c = fh.read(chunksize)
                        if not c:
                            break
                        buf.write(c)
                        buf.seek(0, 0)
                        st = obspy.read(buf,format='MSEED', starttime=starttime,endtime=endtime)
                    # Do something useful!
                    #print(st)
                    if len(st) > 0:
                        if madeStream:
                            stream=stream+st
                        else:
                            stream=st
                            madeStream=True
                        stream.merge()
                if verbose:
                    print("="*50)
                    print(stream)
        if firstTime:
            outstream=stream
        else:
            outstream=outstream+stream
        firstTime=False
    return outstream