""""
Dirac Comb routines to detect and remove transients

Versions:
"""
import math as M
import time

from obspy.core import Trace
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import convolve
from scipy.fftpack import ifft

from .utils import stack_data
from ..time_spans import TimeSpans
from ..logger import init_logger

logger = init_logger()
DEBUG = False


def comb_calc(inp, tp, plots, noise_spans, slice_starttime):
    """
    Calculate dirac comb and average transient

    Args:
        inp (~class `obspy.core.stream.Trace`): the data to be cleaned
        tp: Transient parameters
        plots (bool): produce plots
        noise_spans (~class .TimeSpans):
        slice_starttime (~class UTCDateTime): start of first transient "slice"

    Returns:
        (tuple)
            (obspy.Trace):averaged transient
            (obspy.Trace): dirac comb
            (int): number of teeth in comb,
            (np.array): transient starting one sample earlier
            (np.array): transient starting one sample later,
            (?): )dirac comb buffer
    """
    print(f"Running comb_calc({tp}")
    # c1, c2 = tp.clips[0], tp.clips[1]
    sps = inp.stats.sampling_rate
    dt = 1 / sps
    noise_spans.append(
        _remove_noisy(inp, noise_spans, slice_starttime, tp.period, plot=plots)
    )
    # eq_template = _remove_noisy(inp, eq_template, slice_starttime,
    #                             tp.period, plot=plots)
    tp.period = _refine_period(tp, inp, noise_spans, slice_starttime)

    # Calculate signal minus transients
    x, comb, cbuff, nt, xm, xp, xmm, xpp = comb_stack(
        inp, tp.period, plots, noise_spans, slice_starttime, tp.clips
    )
    # print(f'comb_calc(): {cbuff=}')

    # Warn if the transient comes too late in the period (should not happen)
    imax = x.argmax()
    if imax * dt > tp.period / 2:
        print(
            "Transient may be too close to end of period (peak at "
            "{:.0f}% of test_per)".format(100 * imax * dt / tp.period)
        )
        print("check figure 103 for position\n.")
    dirac_comb = inp.copy()
    dirac_comb.data = comb
    transient = Trace(x)
    transient.stats.sampling_rate = inp.stats.sampling_rate
    return transient, dirac_comb, nt, xm, xp, cbuff


def comb_remove(inp, tp, match, slice_starttime, plots=False):
    """
    Remove near-periodic transients from signal

    :param inp: the data to be cleaned
    :type inp: ~class `obspy.core.stream.Trace`
    :param tp: PeriodicTransient object
    :param match: match and cancel each pulse separately
    :type match: bool
    :param slice_starttime: first slice starttime
    :param plots: plot stuff
    :type plots: bool
    :returns: cleaned signal, transient time series used to clean
    :rtype: obspy.trace, obspy.trace
    """
    print(f"Running comb_clean({tp}, match={match}")
    # Calculate corrected data
    assert isinstance(inp, Trace)
    out, synth = _comb_remove_all(
        inp, tp.dirac_comb, tp.transient_model, tp.comb_buffer, plots
    )
    if match:
        out = _match_each(inp, out, synth, slice_starttime, tp, plots=plots)

    return out, synth


def _refine_period(tp, inp, noise_spans, slice_starttime):
    """
    Find the best average period between transients (if tp.dp > 0)

    Args:
        tp (~class `.TransientParameter`):
        inp (~class `obspy.core.stream.Trace`): input data
        noise_spans (~class `TimeSpans`):
        slice_starttime (~class `UTCDateTime`): start time of first slice
    """
    if not tp.dp > 0:
        return tp.period
    else:
        yy = np.zeros((2, 3))
        c1, c2 = tp.clips[0], tp.clips[1]
        for i in range(3):
            per = tp.period + tp.dp * (i - 1)
            yy[0, i] = per
            x, c, cbuff, *_ = comb_stack(
                inp, per, False, noise_spans, slice_starttime, tp.clips
            )
            z, _ = _comb_remove_all(inp, c, x, cbuff)
            # z, *_ = comb_stack(inp, tp, False, noise_spans, slice_starttime,
            #                    tp.clips)
            # print(noise_spans)
            yy[1, i] = np.sum(
                noise_spans.zero(z).data.clip(c1, c2) ** 2
            ) / np.sum(noise_spans.zero(inp).data.clip(c1, c2) ** 2)
        eps = (yy[1, 0] - yy[1, 2]) / (yy[1, 0] + yy[1, 2] - 2 * yy[1, 1]) / 2
        tp.period += eps * tp.dp
        print("\tbest period found={:g}".format(tp.period))
        if abs(eps) > 1:
            print(">>> outside of test_period+-dp: increase dp2 ", end="")
            print("or change your interval to {:g} <<<".format(tp.period))
        return tp.period


def _match_each(
    inp, out, synth, slice_starttime, tp, adjust_limit=2, plots=False
):
    """
    Individually shift each transient to best match data

    :param inp: input waveform
    :param out: output waveform
    :param synth: synthetic transient waveform
    :param slice_starttime: starttime of first slice
    :param tp: PeriodicTransient object

    SHOULD WE CALCULATE c USING CLIPPED DATA?
       PRO: KEEPS MATCH CORRECTION VALUES FROM GOING CRAZY IN EARTHQUAKES
       CON: MAKES SURE TRANSIENT IS IN DATA (EARTHQUAKES MAY PUMP TRANSIENT
            BEYOND CLIPPED RANGE
    ALTERNATIVE: IF MATCH CORRECTION VALUES GO BEYOND -1,1 DON'T USE THEM?
    """
    dt = 1 / out.stats.sampling_rate
    start_addr = M.floor((slice_starttime - inp.stats.starttime) / dt)
    assert start_addr >= 0
    # nx = xm.size

    print("\tIndividually matching transients")
    k = M.floor((out.stats.npts - start_addr) / (tp.period / dt))
    # data_clipped = inp.data.clip(tp.clips[0], tp.clips[1])
    for i in range(k):  # 0,1...k-1   =1:k
        out = _match_one(
            out, i, synth, slice_starttime, tp, adjust_limit, plots
        )

    # I COULD (SHOULD?) RECALCULATE THE TRANSIENT USING THE IMPROVED COMB
    print("Done individually matching transients")
    if plots:
        hours = np.arange(inp.stats.npts) / (3600 / dt)
        plt.figure(105)
        plt.plot(hours, inp.data, "b", label="signal")
        plt.plot(hours, synth.data, "r", label="synthetic")
        plt.legend()
        plt.xlabel("hours")
        plt.title("Signal and matched synthetic transient")
        plt.show()
    return out


def _match_one(out, i, synth, slice_starttime, tp, adjust_limit=2, plot=False):
    """
    Find the best combination of the transient shifted one left and one right

    Find [amp_left;amp_right] that minimize energy in the equation
        [transient_left(:) transient_right(:)]*[amp_left; amp_right] = out

    :param out: data trace to be corrected
    :param synth: synthetic transient waveform
    :param slice_starttime: starttime of first slice
    :param tp: PeriodicTransient object
    """
    dt = 1 / out.stats.sampling_rate
    nx = tp.transient_model.data.size
    # print(f'{nx=}')
    n = out.stats.npts
    # round to nearest integer
    n1 = M.floor((slice_starttime - out.stats.starttime + i * tp.period) / dt)
    n2 = n1 + nx
    assert n1 > 0
    assert n2 <= n, f"{n2=} is beyond end of data ({n=})"
    # Find best combination of one before and one after to match transient
    g = np.zeros((nx, 1))
    g[:, 0] = out.data[n1:n2]  # Clipped version for calc
    A = np.array([tp.tm, tp.tp]).T
    # Solve c for A*c = g
    c, res, rank, s = np.linalg.lstsq(A, g, rcond=None)
    # IGNORE CORRECTIONS THAT ARE TOO LARGE
    if np.mean(abs(c)) > adjust_limit:
        print(f"Did not tune transient {i:d}: mean of adjustments > allowed:")
        print(
            "\tmean(abs([{:.3g}, {:.3g}])) > {:g}".format(
                c[0, 0], c[1, 0], adjust_limit
            )
        )
    res_transient = np.matmul(A, c)
    gg = g - res_transient

    # The [:,0]s at the end reduce the arrays from Mx1 to one-dim
    # M-length
    comb = tp.dirac_comb.data
    out.data[n1:n2] = gg[:, 0]
    synth.data[n1:n2] = synth.data[n1:n2] + res_transient[:, 0]
    # Set up comb shifted one to right ("m") and one to left ("p")
    # This seems contradictory, but using [n1-1:n2-1] ("minus") shifts
    # delta to the right and is consistent with definitions of tp.tm & tp.tp
    if n1 == 0:
        comb_m = np.hstack((np.zeros(1), comb[n1: n2 - 1]))
    else:
        comb_m = comb[n1 - 1: n2 - 1]
    if n2 == n:
        comb_p = np.hstack((comb[n1 + 1: n2], np.zeros(1)))
    else:
        comb_p = comb[n1 + 1: n2 + 1]
    comb_shifts = np.array([comb_m, comb_p]).T
    shifted_comb = comb[n1: n2] + np.matmul(comb_shifts, c)[:, 0]
    if plot:
        plt.figure(51)
        nn = 10  # number of samples to plot on each side of peak
        iMax = comb[n1: n2].argmax()
        plt.plot(
            comb[n1 + iMax - nn: n1 + iMax + nn],
            "r",
            linewidth=3,
            label="original",
        )
        plt.plot(
            comb_m[iMax - nn: iMax + nn] * c[0, 0], "r--", label="modifiers"
        )
        plt.plot(comb_p[iMax - nn: iMax + nn] * c[1, 0], "r")
        plt.plot(
            shifted_comb[iMax - nn: iMax + nn],
            "b",
            linewidth=3,
            label="matched",
        )
        plt.legend()
        title_text = "{:d}:[[{:.3f}] [{:.3f}]]".format(i, c[1, 0], c[0, 0])
        plt.title(title_text)
        print(title_text)
        plt.show()
    return out


def comb_stack(inp, period, plots, noise_spans, first_slice_starttime, clip):
    """
    Calculate periodic transients and remove from signal

    Args:
        inp (~class `obspy.core.stream.Trace`): the series to be cleaned
        period (float): the transient period (in seconds)
        plots (bool): produce plots
        noise_spans (~class `.TimeSpans`):
        first_slice_starttime (~class `obspy.core.UTCDateTime`): start the
            dirac comb at this time
        clip (tuple): low and high values to clip trace at when calculating
            average transient

    Returns:
        (tuple):
            x   (np.array)   :  averaged transient
            c  (np.array)    :  comb
            cbuff: comb buffer
            ng  (float)      :  number of transients used to make average
            xm  (np.array)   :  x starting one sample earlier
            xp  (np.array)   :  x starting one sample later
            xmm  (np.array)   :  x starting two samples earlier
            xpp  (np.array)   :  x starting two samples later

    xm, xp, xmm and xpp should probably be removed and replaced by in-place
    shifting of x (with simple padding at the ends)
    """
    print(f"\tRunning comb_stack({period=}s)", flush=True)
    dt = inp.stats.delta
    samps_per_period = period / dt
    n = inp.stats.npts
    rp = int(np.round(samps_per_period))

    # -----------------------------------------
    # Create the Dirac comb, starting at slice_starttime
    slice_startaddr = M.floor(
        (first_slice_starttime - inp.stats.starttime) / dt
    )
    assert slice_startaddr >= 0
    if DEBUG:
        print("COMB_STACK(): Creating Dirac comb")
    n_teeth, cbuff, c = comb(n - slice_startaddr, samps_per_period)
    # print(f'comb_stack(): {cbuff=}')
    if not slice_startaddr == 0:
        # Stuff slice_startaddr zeros before dirac
        c = np.hstack((np.zeros(slice_startaddr), c))

    # Create the broken Dirac comb (teeth missing for slices containing
    if DEBUG:
        print(
            "COMB_STACK(): Creating broken Dirac comb...", flush=True, end=""
        )
        tic = time.time()
    cb = c.copy()
    b_teeth = n_teeth
    # if plots:
    #    eq_template.plot_mask(inp)
    for i in range(n_teeth):
        slice_starttime = first_slice_starttime + i * period
        slice_endtime = slice_starttime + period
        n1 = int(slice_startaddr + M.floor(i * samps_per_period))
        n2 = n1 + rp
        # if np.any(eq_template.data[n1: n2] == 0):
        if noise_spans.has_zeros(slice_starttime, slice_endtime):
            cb[n1:n2] = 0
            b_teeth -= 1
    if DEBUG:
        print("Took {:.1f} s".format(time.time() - tic))

    if plots:
        plt.figure(102)
        plt.subplot(2, 1, 1)
        plt.plot(np.arange(c.size) * dt, c)
        plt.title("comb")
        plt.subplot(2, 1, 2)
        plt.plot(np.arange(cb.size) * dt, cb)
        plt.title("broken comb")
        plt.xlabel("Time (seconds)")
        plt.show()
    x, xm, xp, xmm, xpp = _calc_transient(
        inp, cbuff, rp, cb, b_teeth, slice_startaddr, clip, plots
    )
    return x, c, cbuff, b_teeth, xm, xp, xmm, xpp


def _comb_remove_all(
    inp, dirac_comb, transient_model, comb_buffer, plots=False
):
    """
    Remove periodic transients from signal

    Args:
        inp (obspy.core.Trace): the Trace to be cleaned
        dirac_comb (obspy.core.Trace): dirac comb (one spike for each
            transient)
        transient_model (obspy.core.Trace): averaged transient
        comb_buffer (int): number of buffer samples inserted before true
            comb start
        tp (Periodic_Transient): holds dirac_comb and transient_model

    Returns:
        (tuple):
            (obspy.Trace): cleaned series (inp - y)
            (obspy.Trce): sythetic transient series (y)
    """
    print("\tRunning _comb_remove_all", flush=True)
    n = inp.stats.npts

    # Create synthetic transient time series
    if DEBUG:
        print("COMB_STACK(): Creating synthetic transient time series")
    y = inp.copy()
    # print(f'{comb_buffer=}')
    yy = convolve(dirac_comb.data, transient_model.data)
    y.data = yy[comb_buffer: comb_buffer + n]

    if plots:
        dt = inp.stats.delta
        index = np.arange(n) / dt
        plt.figure("_comb_remove_all()")
        plt.plot(index, inp.data, "b", label="signal")
        plt.plot(index, y.data, "r", label="synthetic")
        plt.legend()
        plt.xlabel("seconds")
        plt.title("Signal and synthetic transient")
        plt.show()

    # Create corrected time series
    if DEBUG:
        print("COMB_STACK(): Creating corrected time series")
    out = inp.copy()
    out.data = inp.data - y.data

    return out, y


def _calc_transient(
    inp, cbuff, rp, cb, b_teeth, slice_startaddr, clip, plots=False
):
    """
    Calculate the average transient by correlating around the broken comb

    :param inp: input data (Trace)
    :param cbuff: comb buffer (calculated by comb())
    :param rp: number of samples on either side of comb tooth to use
    :param cb: broken Dirac comb
    :param b_teeth: teeth in broken comb
    :param clip: clip levels (2-tuple)
    :param slice_startaddr: address of first slice start
    """
    # ** Could make this faster by reducing the buffer on either
    # ** side to that which is extracted in the next section

    # adding x zeros to either side of the comb gives a length 2*rp+1 result
    # from correlate('valid')
    # buff=np.hstack( (np.zeros(rp),cb,np.zeros(rp)) )
    # Instead of buffering by rp (above), why not just buffer by a little more
    # than cbuff (all we recover is ~-cbuff-2:rp-cbuff+2)
    nbuff = cbuff + 5
    buff = np.hstack((np.zeros(rp), cb, np.zeros(nbuff)))
    if DEBUG:
        print("_calc_transient(): Calculating average transient by"
              "correlation...", flush=True, end="")
        tic = time.time()
    xx = (
        np.correlate(inp.data.clip(clip[0], clip[1]), buff, mode="valid")
        / b_teeth
    )
    if DEBUG:
        print("Took {:.1f} s".format(time.time() - tic))

    # Extract average transient and same shifted one sample right, left
    if DEBUG:
        print(f"{xx.size=}, {rp=}, {cbuff=}")
        print("_calc_transient(): Extracting average transient")
    x = xx[nbuff - cbuff: nbuff + rp - cbuff + 1]
    xm = xx[nbuff - cbuff - 1: nbuff + rp - cbuff]
    xp = xx[nbuff - cbuff + 1: nbuff + rp - cbuff + 2]
    xmm = xx[nbuff - cbuff - 2: nbuff + rp - cbuff - 1]
    xpp = xx[nbuff - cbuff + 2: nbuff + rp - cbuff + 3]
    if plots:
        dt = 1.0 / inp.stats.sampling_rate
        hours = np.arange(rp) / (3600 / dt)
        plt.figure(num="calc_transient")
        plt.plot(
            hours,
            inp.data[slice_startaddr: slice_startaddr + rp],
            "b",
            label="Example signal",
        )
        plt.plot(hours, x[:rp], "r", label="comb_stack", linewidth=2)
        plt.xlabel("hours")
        plt.title("blue: signal, red: comb_stack")
        plt.legend()
        plt.show()
    return x, xm, xp, xmm, xpp


def comb(n, per):
    """
    Generate a Dirac comb with n/per "teeth" (impulses)

    Set up in freq domain, then calculate the ifft

    WOULD PROBABLY BE MUCH FASTER IF CALCULATED FOR pow2(N), then cut down the
    RESULT TO LENGTH N

    :param n: length (samples)
    :param per: period (in samples)
    Output:
        n_teeth - number of teeth created
        cbuff    - offset before teeth
        out    - Dirac comb
    """
    # print('Running comb(n={:g},per={:g})'.format(n,per),flush=True)
    if DEBUG:
        print("COMB(): Generating DIRAC comb")

    # CHANGE N TO A POWER OF TWO
    n_return = n
    n = 2 ** (M.ceil(M.log2(n_return)))
    if DEBUG:
        print("COMB(): n={:d}, npow2={:d}".format(n_return, n))

    n_teeth = M.ceil(n / per)
    cl = np.zeros(n) + 1j * np.zeros(n)
    # teeth are shifted "cbuff" samples to the right (to eliminate edge
    # effects?) to be removed after the convolution in comb_stack
    cbuff = 48
    offex = 2 * M.pi * 1j * cbuff / n
    fn2 = M.floor(n / 2)
    fn4 = M.floor(n / 4)
    if DEBUG:
        print(
            "COMB(): Filling positive freq terms "
            "(0:{:d}) with fft of deltas...".format(fn2),
            flush=True,
            end="",
        )
        tic = time.time()

    # ORIGINAL METHOD
    #     for l in range(fn2+1) :  # for l=0:fn2
    #         lpn=l*per/n;
    #         if np.abs(np.round(lpn)-lpn)<1e-9 :
    #             cl[l]=n_teeth;
    #         else :
    #             q=np.exp(-2*M.pi*1j*lpn)
    #             cl[l]=(q**n_teeth-1)/(q-1)
    #         cl[l]=cl[l]*np.exp(-offex*l)
    # ALTERNATIVE (FASTER) METHOD
    lf = np.arange(fn2 + 1)
    lpn = lf * per / n
    q = np.exp(-2 * M.pi * 1j * lpn)
    iPulses = np.abs(np.round(lpn) - lpn) < 1e-9
    q[iPulses] = -1
    cl[: fn2 + 1] = np.exp(-offex * lf) * (q**n_teeth - 1) / (q - 1)
    cl[: fn2 + 1][iPulses] = n_teeth
    if DEBUG:
        print("Took {:.1f} seconds".format(time.time() - tic))
        print(
            "COMB(): multiplying upper half of frequencies "
            "({:d}:{:d}) by cos(0->pi)...".format(fn4, fn2),
            flush=True,
            end="",
        )
        tic = time.time()
    for ff in range(fn4, fn2 + 1):  # l=fn4:fn2
        cl[ff] *= (1 + M.cos((ff - fn4) / (fn2 - fn4) * M.pi)) / 2.0
    if fn2 == n / 2:  # even # of elements
        # cl(fn2+2:n)=conj(cl(fn2:-1:2));
        cl[fn2 + 1: n] = np.conj(cl[fn2 - 1: 0: -1])
    else:
        # cl(fn2+2:n)=conj(cl(fn2+1:-1:2));
        cl[fn2 + 1: n] = np.conj(cl[fn2: 0: -1])
    if DEBUG:
        print("Took {:.1f} seconds".format(time.time() - tic))
        print(
            "COMB(): Calculating IFFT of {:d}-point series...".format(cl.size),
            flush=True,
            end="",
        )
        tic = time.time()

    out = np.real(ifft(cl))
    if DEBUG:
        print("Took {:.1f}s".format(time.time() - tic))

    # CONVERT BACK TO ORIGINAL LENGTH
    n_teeth = M.ceil(n_return / per)
    out = out[:n_return]
    return n_teeth, cbuff, out


def _remove_noisy(inp, eq_spans, first_slice_starttime, period, plot=False):
    """
    Returns the time spans of data with anomalously high noise

    Affects indices corresponding to the middle half of a slice (based on
    first_slice_starttime and period), to allow for subsequent small changes
    in period

    Args:
        inp (obspy.trace): the original data
        eq_spans (TimeSpans): EQ time spans to remove
        first_slice_starttime (UTCDateTime): the offset (seconds) of the first
            slice
        period (float): the interval at which the data will be cut into slices
    Returns:
        bad_spans (~class `.TimeSpans`): noisy time spans
    """
    dt = inp.stats.delta
    start_addr = M.floor((first_slice_starttime - inp.stats.starttime) / dt)
    assert start_addr >= 0
    samp_per = period / dt
    # rp = np.round(samp_per)
    n = inp.stats.npts - start_addr
    n_slices = M.floor(n / samp_per)

    # Calculate the variance of each data slice
    slices = list()
    for i in range(n_slices):
        slice_starttime = first_slice_starttime + i * period
        slice_endtime = slice_starttime + period
        # n1 = int(start_addr + M.floor(i*samp_per))
        # n2 = int(n1+rp)
        # Don't go past end
        # if slice_endtime >= inp.stats.npts:
        if slice_endtime >= inp.stats.endtime:
            continue
        # Ignore the first & last 10 samples to accomodate small changes in
        # period
        # var = sum((inp.data[n1+10: n2-10]
        #            * eq_template.data[n1+10: n2-10])**2)
        slice = inp.slice(slice_starttime + 10 * dt, slice_endtime - 10 * dt)
        var = sum(eq_spans.zero(slice).data ** 2)
        hz = eq_spans.has_zeros(slice_starttime, slice_endtime)
        slices.append(
            {
                "var": var,
                "st": slice_starttime,
                "et": slice_endtime,
                "hasZeros": hz,
            }
        )
        # slices.append({'var': var, "n1": n1, "n2": n2, 'hasZeros': hz})

    # Calculate the median and std of variances of slices without zeros
    variances = [x["var"] for x in slices if not x["hasZeros"] is True]
    median = np.median(variances)
    sigma = np.std(variances)

    # For slices with abnormally high variance, set eq_template indices to zero
    print(
        "\tRejecting high-variance slices (>{:.4g}+3*{:.4g})...".format(
            var, sigma
        ),
        end="",
    )
    nRejected = 0
    noise_spans = []
    for slice in slices:
        # print(f'{slice["var"]=}')
        if slice["var"] <= median + 3 * sigma:
            slice["muted"] = False
        else:
            nRejected += 1
            slice["muted"] = True
            noise_spans.append([slice["st"] + 10 / dt,
                                slice["et"] - 10 / dt])
            # eq_template.data[n1+10:n2-10] = 0
    if nRejected > 0:
        print("{:d} of {:d} rejected".format(nRejected, len(slices)))
    else:
        print("none rejected")

    if plot:
        stack = stack_data(inp.data[start_addr:], samp_per)
        _plot_remove_noisy(stack, slices, dt)
    return TimeSpans(noise_spans)


def _plot_remove_noisy(stack, slices, dt):
    plt.figure(num="remove_noisy")
    nrows, ncols = stack.shape
    # print(f'{nrows=}, {ncols=}, {len(slices)=}')
    time = np.arange(nrows) * dt
    hasZeros = [i for i in range(ncols) if slices[i]["hasZeros"] is True]
    isMuted = [i for i in range(ncols) if slices[i]["muted"] is True]
    isUsed = [
        i
        for i in range(ncols)
        if not slices[i]["hasZeros"] is True and not slices[i]["muted"] is True
    ]
    # print(hasZeros)
    # print(isMuted)
    lines, labels = [], []
    if len(hasZeros) > 0:
        ls = plt.plot(time, stack[:, hasZeros], "r:")
        lines.append(ls[0])
        labels.append("Has EQ Zeros")
    if len(isMuted) > 0:
        ls = plt.plot(time, stack[:, isMuted], "g--")
        lines.append(ls[0])
        labels.append("High Variance")
    if len(isUsed) > 0:
        ls = plt.plot(time, stack[:, isUsed], "b-")
        lines.append(ls[0])
        labels.append("Clean (used)")
    plt.legend(lines, labels)
    plt.title(
        "{} Slices: {} EQZeroed, {} high Variance (muted), {} Used".format(
            len(slices), len(hasZeros), len(isMuted), len(isUsed)
        )
    )
    plt.show()
