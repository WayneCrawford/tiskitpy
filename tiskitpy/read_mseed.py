#!env python3
""" Read miniSEED files, handling files that are too big for obspy"""
import obspy
import glob
import io
import os.path

from tiskitpy.logger import init_logger

logger = init_logger()

def read_MSEED(filenames, starttime, endtime, verbose=False):
    """
    Read in miniSEED files, using either obspy.read or, if the file is
    more than 2.1 Gb (libmseed limitation), by reading in chunks first
    """
    firstTime = True
    for filename in glob.glob(filenames):
        filesize = os.path.getsize(filename)
        if filesize < 2**31:
            if verbose:
                print("Filesize {:.2f}Go < 2^31 bytes, reading using "
                      "obspy.read".format(filesize/float(1024*1024*1024)))
            stream = obspy.read(filename, format='MSEED',
                                starttime=starttime,
                                endtime=endtime)
        else:
            if verbose:
                print("Filesize {:.2f}Go >= 2^31 bytes, breaking into "
                      "readable chunks for obspy.read"
                      .format(filesize/float(1024*1024*1024)))
            reclen = 4096
            chunksize = 20000 * reclen  # Around 80 MB
            madeStream = False

            with io.open(filename, "rb") as fh:
                while True:
                    with io.BytesIO() as buf:
                        c = fh.read(chunksize)
                        if not c:
                            break
                        buf.write(c)
                        buf.seek(0, 0)
                        st = obspy.read(buf, format='MSEED',
                                        starttime=starttime,
                                        endtime=endtime)
                    # Do something useful!
                    if len(st) > 0:
                        if madeStream:
                            stream = stream + st
                        else:
                            stream = st
                            madeStream = True
                        stream.merge()
                if verbose:
                    print("="*50)
                    print(stream)
        if firstTime:
            outstream = stream
        else:
            outstream = outstream + stream
        firstTime = False
    return outstream
