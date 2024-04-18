#!/usr/bin/python3
"""
Convert an obsPy data stream from zero phase to minimum phase

Duplicates Scherbaum C programs "interpolate", "fir2caus" and "decimate"

Uses a precalculated FIR "correction" file made using the Matlab code
makeFIRcorrFile.m (or imported from Scherbaum FIR_CORR directory) and then
converted to JSON using my conv_ScherbaumParm2JSON.py code

Author: Wayne Crawford
"""

import numpy as np
import json as json
import scipy.signal as sig

from ..logger import init_logger

logger = init_logger()



def fir2caus(stream, firCorrFileName, FIRdecim):
    """
    Main conversion routine, calls all others

    Args:
        stream (:class:`obspy.core.stream.Stream`): input waveforms
        firCorrFileName (str): name of JSON file containing FIR corrections
        FIRdecim (int): decimation factor associated with FIR to correct
    Returns:
        s (:class:`obspy.core.stream.Stream`): output waveforms
    """
    ##################################################

    s = stream.copy()
    for tr in s:
        data_type = tr.data.dtype
        sample_frequency_out = tr.stats.sampling_rate
        # Extract data  (1D numpy array)
        tr.detrend()
        data = tr.data

        # INTERPOLATE
        data_interp = _interpolate(data, FIRdecim, False)
        # CHANGE FIR FROM CAUSAL TO MIN PHASE
        data_interp_corr, offset_npts = fir2caus_subr(
            data_interp, firCorrFileName, False
        )
        # RESAMPLE TO ORIGINAL
        data_corr = data_interp_corr[::FIRdecim]

        tr.data = np.array(
            data_corr, dtype=data_type
        )  # Force back to original data type:
        tr.stats.starttime += offset_npts / (sample_frequency_out * FIRdecim)
    return s


def _interpolate(trace, ipol_fac, returnOrigDType=True):
    """
    Interpolate data by an integer factor that is a multiple of 2,3, and/or 5

    Args:
        trace (:class:`numpy.ndarray`): data to interpolate
        ipol_fac (int): factor to interpolate by
        returnOrigDType (bool): return output trace as same type as input
            trace (True) or as float double (False)
    """

    if ipol_fac == 1:
        print("Nothing to interpolate, exit...\n")
        return trace

    conv_fac = int(ipol_fac)

    # number of divisions by 2
    pow2 = 2
    no_it2 = 0
    while conv_fac / pow2 == int(conv_fac / pow2):
        pow2 *= 2
        no_it2 += 1
    if no_it2 > 0:
        conv_fac = 2 * conv_fac / pow2

    # number of divisions by 3 */
    pow3 = 3
    no_it3 = 0
    while conv_fac / pow3 == int(conv_fac / pow3):
        pow3 *= 3
        no_it3 += 1
    if no_it3 > 0:
        conv_fac = 3 * conv_fac / pow3

    #  number of divisions by 5 */
    pow5 = 5
    no_it5 = 0
    while conv_fac / pow5 == int(conv_fac / pow5):
        pow5 *= 5
        no_it5 += 1
    if no_it5 > 0:
        conv_fac = 5 * conv_fac / pow5

    if conv_fac != 1:
        print(
            "Factor {:d} cannot be separated in factors 2, 3, and 5\n".format(
                ipol_fac
            )
        )
        exit(1)

    data_type = trace.dtype
    # interpolations by factors of  2
    for k in range(0, no_it2):  # (k=0;k<no_it2;k++) :
        trace = _ipol2(trace)
    # interpolations by factors of 3
    for k in range(0, no_it3):  # (k=0;k<no_it3;k++) :
        trace = _ipol3(trace)
    # interpolations by factors of  5
    for k in range(0, no_it5):  # (k=0;k<no_it5;k++)
        trace = _ipol5(trace)

    if returnOrigDType:
        trace = np.array(
            trace, dtype=data_type
        )  # Force back to original data type:
    else:
        trace = np.array(trace)

    return trace


# E. Wielandt's filter coefficients for interpolation ratios 2, 3 and 5 */
# g2_1 RMS rel. error for k_rel <0.8: 0.009% */
# g3_1 RMS rel. error for k_rel <0.8: 0.008% */
# g5_1 RMS rel. error for k_rel <0.8: 0.006% */
# g5_2 RMS rel. error for k_rel <0.8: 0.009% */


def _ipol2(xin, debug=False):
    "interpolation by factor 2"

    g2_1 = np.array(
        [
            -0.0002616,
            0.0009302,
            -0.0023258,
            0.0049142,
            -0.0092537,
            0.0161355,
            -0.0266150,
            0.0424554,
            -0.0670612,
            0.1091686,
            -0.2008408,
            0.6327541,
            0.6327542,
            -0.2008408,
            0.1091686,
            -0.0670612,
            0.0424554,
            -0.0266150,
            0.0161355,
            -0.0092537,
            0.0049142,
            -0.0023258,
            0.0009302,
            -0.0002616,
        ],
        dtype="double",
    )

    # xx = interpolated sample at center between index i and i+1
    xx = np.convolve(xin, g2_1[::-1], mode="full")
    if debug:
        print(f"len(g2_1)={len(g2_1)}")
    xx = xx[int(len(g2_1) / 2): (len(xin) + int(len(g2_1) / 2))]
    # print("{:d} {:d}".format(len(xin),len(xx)))
    # WEAVE data together
    xout = np.column_stack((xin, xx))  # 2D
    xout = xout.flatten()  # To 1D
    xout = xout[0:-1]  # get rid of extrapolations
    # Should verify that xout[0]==xin[0] and that xout[-1]==xin[-1]
    return xout


def _ipol3(xin):
    """interpolation by factor 3"""
    ifac = 3  # interpolation factor
    g3_1 = np.array(
        [
            -0.0002324,
            0.0008256,
            -0.0020654,
            0.0043684,
            -0.0082410,
            0.0144087,
            -0.0238643,
            0.0383080,
            -0.0611438,
            0.1015046,
            -0.1959616,
            0.8226883,
            0.4111037,
            -0.1564917,
            0.0885536,
            -0.0553522,
            0.0353777,
            -0.0223077,
            0.0135759,
            -0.0078056,
            0.0041522,
            -0.0019670,
            0.0007871,
            -0.0002211,
        ],
        dtype="double",
    )

    # xx1: interpolated sample 1/3 between index i & i+1
    # xx2: interpolated sample 2/3 between index i & i+1 : use g3_1[] backwards
    xx1 = np.convolve(xin, g3_1[::-1], mode="full")
    xx2 = np.convolve(xin, g3_1, mode="full")
    xx1 = xx1[len(g3_1) / 2: (len(xin) + len(g3_1) / 2)]
    xx2 = xx2[len(g3_1) / 2: (len(xin) + len(g3_1) / 2)]
    # WEAVE xx1 and xx2 into xin
    xout = np.column_stack((xin, xx1, xx2))
    xout = xout.flatten()
    xout = xout[0: (-ifac + 1)]  # get rid of extrapolations
    # Should verify that xout[0]==xin[0] and that xout[-1]==xin[-1]
    return xout


def _ipol5(xin):
    """interpolation by factor 5"""
    ifac = 5  # interpolation factor
    g5_1 = np.array(
        [
            -0.0001611,
            0.0005720,
            -0.0014316,
            0.0030307,
            -0.0057263,
            0.0100352,
            -0.0166788,
            0.0269186,
            -0.0433567,
            0.0732492,
            -0.1480766,
            0.9320452,
            0.2327664,
            -0.0984035,
            0.0572469,
            -0.0362360,
            0.0233232,
            -0.0147709,
            0.0090151,
            -0.0051932,
            0.0027660,
            -0.0013112,
            0.0005248,
            -0.0001473,
        ],
        dtype="double",
    )
    g5_2 = np.array(
        [
            -0.0002526,
            0.0008977,
            -0.0022452,
            0.0047467,
            -0.0089479,
            0.0156272,
            -0.0258392,
            0.0413715,
            -0.0657512,
            0.1082703,
            -0.2048060,
            0.7525200,
            0.5015040,
            -0.1790148,
            0.0997642,
            -0.0619420,
            0.0394431,
            -0.0248145,
            0.0150789,
            -0.0086611,
            0.0046043,
            -0.0021804,
            0.0008723,
            -0.0002452,
        ],
        dtype="double",
    )

    # xx1 = interpolated sample 1/5 between index i & i+1 : use g5_1
    # xx2 = interpolated sample 2/5 between index i & i+1 : use g5_2
    # xx3 = interpolated sample 2/5 between index i & i+1 : use g5_2 backwards
    # xx4 = interpolated sample 2/5 between index i & i+1 : use g5_1 backwards
    xx1 = np.convolve(xin, g5_1[::-1], mode="full")
    xx2 = np.convolve(xin, g5_2[::-1], mode="full")
    xx3 = np.convolve(xin, g5_2, mode="full")
    xx4 = np.convolve(xin, g5_1, mode="full")
    xx1 = xx1[len(g5_1) / 2: (len(xin) + len(g5_1) / 2)]
    xx2 = xx2[len(g5_2) / 2: (len(xin) + len(g5_2) / 2)]
    xx3 = xx3[len(g5_2) / 2: (len(xin) + len(g5_2) / 2)]
    xx4 = xx4[len(g5_1) / 2: (len(xin) + len(g5_1) / 2)]

    # WEAVE data together
    xout = np.column_stack((xin, xx1, xx2, xx3, xx4))
    xout = xout.flatten()
    xout = xout[0: (-ifac + 1)]  # get rid of extrapolations
    # Should verify that xout[0]==xin[0] and that xout[-1]==xin[-1]
    return xout


def fir2caus_subr(intrace, firCorrFileName, returnOrigDType=True):
    """
    Implements FIR filter correction

    Described in chapter 8 of Scherbaum, 1996, "Of poles and zeros"

    Args:
        intrace (:class:`numpy.ndarray`: 1d float array
        firCorrFileName (str): JSON parameter file name
        returnOrigDtype (bool): return output trace as same type as input
            trace (True) or as float double (False)

    Returns:
        (tuple):
            outtrace (:class:`numpy.ndarray`): causal version of trace
            timetag (int): number of samples that outtrace delta response is
                BEFORE intrace delta response


    BASED ON Scherbaum C code:
    VERSION: 2.0
    DATE: 1996-10-18 (Frank Scherbaum)
    DESCRIPTION: This program implements the FIR filter correction described in
    chapter 8 of Scherbaum, F., Of poles and zeros: fundamentals of digital
    seismology, Kluwer Academic Publishers, 1996. The key equation
    is (8.15) on page 120. The ARMA coefficients are assumed to be
    stored in a '*.prt' file produced by Frank Scherbaum's analyse.m
    Mathematica program. This is the  <filter coefficient file
    for correction filter> above.

    Filter x[] which is the reversed input sequence x2[]
    using the difference equation:

            mx                mx
            --                -----
    y'[i] =  > a[k]*y'[i-k]   + > b[l] x[i-k]
            __                __
            k=1               l=0

    This corresponds to equ. (8.15) in 'Scherbaum, F: Of poles and zeros,
    Fundamentals of Digital Seismology, Kluwer Academic Publ., 1996'
    mx = number of AR coefficients
    b[l] = MA coefficients for l = 0, mx
    a[k] = AR coefficients for k = 1, mx

    x[i] = reversed input sequence for i = 0 .....
    y'[] = output sequence

    Reverse the output sequence y'[] in time again to obtain the
    corrected sequence y[n]!
    """

    data_type = intrace.dtype
    x = np.array(
        intrace, dtype="double"
    )  # Convert to double for manipulations

    # READ CORRECTION FILTER COEEFICIENTS
    with open(firCorrFileName, "r") as f:
        book = json.load(f)
    a = np.array(book["ar"], dtype="double")
    b = np.array(book["ma"], dtype="double")
    # fir_len = len(book["firdata_eff"])
    timetag = float(book["timetag"])
    # fdig = book["fdig"]

    x = x[::-1]  # flip in time
    # ===================
    # Original code
    # ===================
    # mx = len(a)
    # nx = len(x)
    # # Buffer x with mx zeros at beginning
    # x=np.concatenate((np.zeros(mx,dtype='double'),x))
    # # Set up y as all zeros
    # y=np.zeros(len(x),dtype='double')
    # for i in range(0,nx) :                  # (i=0; i< ndat2; i++) {
    #     y[i+mx] = np.dot(x[i+mx:i:-1],b)           # MA
    #     y[i+mx] += np.dot(y[i+mx-1:i:-1],a[1:])    # AR
    # y=y[mx:]  # Remove buffer
    # ===================
    # Faster replacement
    # ===================
    a[1:] *= -1  # because lfilter subtracts "a"s, invert all EXCEPT first
    y = sig.lfilter(b, a, x)

    y = y[::-1]  # flip back in time

    if returnOrigDType:
        yout = np.array(
            y, dtype=data_type
        )  # Force back to original data type:
        return yout, timetag
    else:
        return y, timetag

    return yout, timetag


def main():
    "Main Program"
    print("fir2caus.py is a module, not a runnable program.")
    return 0


if __name__ == "__main__":
    main()
