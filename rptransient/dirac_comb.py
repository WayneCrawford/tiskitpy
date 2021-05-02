""""
Dirac Comb routines to detect and remove glitches

Versions:
  0: straight implementation of Wielandt's Matlab routines
  1: Add features to handle earthquakes and noisy data:
    - Remove "teeth" from the Dirac comb wherever the "template" time series
      covering the same samples as the data contains a zero within the slice
      corresponding to that tooth
    - Also remove teeth for slices that have significantly higher variance than
      the others
    - Use this damaged comb and a clipped version of the data (using the clips
      variable) to generate the average glitch.
    - Apply this average glitch and the full comb to clean the data
  2:
    - In individual matching, allow further shifting of comb tooth:
        - In time
        - In height? (shouldn't allow too much liberty or may create errors)
    - Recalculate glitch based on this comb?
    - Return a sparse dirac comb for the time series, which could be combined
      with others to create a dirac comb for the whole dataset
    - Return the number of teeth used to calculate the master glitch
"""
import sys
import math as M
import time

# import obspy.core
import numpy as np
# import scipy as sp
import matplotlib.pyplot as plt
from scipy.signal import convolve
# from scipy.signal import correlate, deconvolve
from scipy.fftpack import ifft

from .utils import stack_data

DEBUG = False


def comb_clean(inp, tp, match, plots, template, start_offset):
    """
    Remove near-periodic glitches from signal

    :param inp: the data to be deglitched
    :type inp: ~class `obspy.core.stream.Trace`
    :param tp: Transient parameters
    :param match: match and cancel each pulse separately
    :type match: bool
    :param plots: produce plots
    :type plots: bool
    :param template: template of ones and zeros, zero where the
                     comb should not use data for glitch calculation
                     (not used in diracComb_v0)
    :type template: ~class `obspy.core.stream.Trace`
    :param start_offset: seconds to offset the start of the comb so that
                         the glitch will be well placed for modeling
                         (not used in diracComb_v0)

    Outputs:
        out (obspy.trace)        : cleaned signal
        glit_1 (np.ndarray)      : single averaged glitch
        gl_all (obspy.trace)     : glitch time series used to clean
        dirac_comb (obspy.trace) : comb_stack (was sparse in v2)
        nTeeth (int)             : number of teeth in comb (# of glitches avgd)
    """
    print('Running comb_clean(test_period={:g}s,dp={:g}s,match={})'
          .format(tp.period, tp.dp, match), flush=True)
    c1 = tp.clips[0]
    c2 = tp.clips[1]
    sps = inp.stats.sampling_rate
    dt = 1/sps
    template = _remove_noisy(inp, template, start_offset, tp.period)
    # -----------------------------------------
    # Select average period between glitches:
    #   if dp > 0, test neighboring values
    #   otherwise, use provided value
    if tp.dp > 0:
        yy = np.zeros((2, 3))
        for i in range(3):  # [0,1,2]
            tp = tp.period + tp.dp*(i-1)
            yy[0, i] = tp
            z, *_ = comb_stack(inp, tp, False, template, start_offset, tp.clips)
            yy[1, i] = np.sum((z.data*template.data).clip(c1, c2)**2) /\
                np.sum((inp.data*template.data).clip(c1, c2)**2)
        eps = (yy[1, 0]-yy[1, 2]) / (yy[1, 0]+yy[1, 2]-2*yy[1, 1])/2
        tp = tp.period + eps*tp.dp
        print('\tbest period found={:g}'.format(tp))
        if abs(eps) > 1:
            print('>>> outside of test_period+-dp: increase dp2 ', end='')
            print('or change your interval to {:g} <<<'.format(tp))
    else:
        tp = tp.period

    # Calculate signal minus glitches
    out, gl_all, x, xm, xp, comb, nt, xmm, xpp = comb_stack(
        inp, tp, plots, template, start_offset, tp.clips)

    # Warn if the glitch comes too late in the period (should not happen)
    imax = x.argmax()
    if imax*dt > tp.period/2:
        print('Glitch may be too close to end of period (peak at '
              '{:.0f}% of test_per)'.format(100*imax*dt/tp.period))
        print('check figure 103 for position\n.')

    # -----------------------------------------
    # Individually shift each glitch to best match data
    # SHOULD WE CALCULATE c USING CLIPPED DATA?
    #    PRO: KEEPS MATCH CORRECTION VALUES FROM GOING CRAZY IN EARTHQUAKES
    #    CON: MAKES SURE GLITCH IS IN DATA (EARTHQUAKES MAY PUMP GLITCH BEYOND
    #         CLIPPED RANGE
    # ALTERNATIVE: IF MATCH CORRECTION VALUES GO BEYOND -1,1 DON'T USE THEM?
    if match:
        print('\tIndividually matching glitches')
        nx = x.size
        n = out.stats.npts
        k = M.ceil(n/(tp/dt))

        data_clipped = inp.data.clip(tp.clips[0], tp.clips[1])
        # -----------------------------------------
        # For each glitch slice, find the best combination of the glitch
        # shifted one left and one right to further reduce energy
        # by solving for [amp_left;amp_right] in the equation
        #       [glitch_left(:) glitch_right(:)]*[amp_left;amp_right]=out
        for i in range(k):  # 0,1...k-1   =1:k
            adjust_limit = 2.   # Maximum average change to accept
            # round to nearest integer
            n1 = M.floor((i*tp+start_offset) / dt + .5)
            n2 = n1+nx
            # nBuff=M.floor(test_secs*sps)
            if n2 > n+1:
                # Don't try to change beyond the end of the time series
                continue

            # Find best combination of one before and one after to match glitch
            g = np.zeros((nx, 1))
            g[:, 0] = out.data[n1:n2]  # Clipped version for calc
            A = np.array([xm, xp]).T
            # Solve c for A*c = g
            c, res, rank, s = np.linalg.lstsq(A, g)
            # IGNORE CORRECTIONS THAT ARE TOO LARGE
            if np.mean(abs(c)) > adjust_limit:
                print(f'Did not tune glitch {i:d}: mean of', end=' ')
                print('adjustments > allowed:', end=' ')
                print('mean(abs([{:.3g}][{:.3g}])) > {:g}'
                      .format(c[0, 0], c[1, 0], adjust_limit))
                continue
            res_glitch = np.matmul(A, c)
            gg = g-res_glitch

            # The [:,0]s at the end reduce the arrays from Mx1 to one-dim
            # M-length
            out.data[n1:n2] = gg[:, 0]
            gl_all.data[n1: n2] = gl_all.data[n1: n2] + res_glitch[:, 0]
            # Set up comb shifted one to right ("m") and one to left ("p")
            # This seems contradictory, but using [n1-1:n2-1] ("minus") shifts
            # delta to the right and is consistent with definitions of xm & xp
            if n1 == 0:
                comb_m = np.hstack((np.zeros(1), comb[n1:n2-1]))
            else:
                comb_m = comb[n1-1:n2-1]
            if n2 == n:
                comb_p = np.hstack((comb[n1+1: n2], np.zeros(1)))
            else:
                comb_p = comb[n1+1: n2+1]
            comb_shifts = np.array([comb_m, comb_p]).T
            shifted_comb = comb[n1:n2] + np.matmul(comb_shifts, c)[:, 0]
            if False:
                plt.figure(51)
                nn = 10   # number of samples to plot on each side of peak
                iMax = comb[n1: n2].argmax()
                plt.plot(comb[n1+iMax-nn: n1+iMax+nn], 'r', linewidth=3,
                         label="original")
                plt.plot(comb_m[iMax-nn: iMax+nn]*c[0, 0], 'r--',
                         label="modifiers")
                plt.plot(comb_p[iMax-nn: iMax+nn]*c[1, 0], 'r')
                plt.plot(shifted_comb[iMax-nn:iMax+nn], 'b', linewidth=3,
                         label="matched")
                plt.legend()
                title_text = '{:d}:[[{:.3f}] [{:.3f}]]'.format(
                    i, c[1, 0], c[0, 0])
                plt.title(title_text)
                print(title_text)
                plt.show()
        # end for i in range(k)

        # I COULD (SHOULD?) RECALCULATE THE GLITCH USING THE IMPROVED COMB

        print('Done individually matching glitches')
        if plots:
            hours = np.arange(n) / (3600/dt)
            plt.figure(105)
            plt.plot(hours, inp.data, 'b', label='signal')
            plt.plot(hours, gl_all.data, 'r', label='synglitch')
            plt.legend()
            plt.xlabel('hours')
            plt.title('Signal and matched synthetic glitch')
            plt.show()
    # end "if match :"
    # --------------------------------------------
    # Significant change in v3: use obspy trace instead of sparce np.array
    dirac_comb = out.copy()
    dirac_comb.data = comb

    return (out, x, gl_all, dirac_comb, nt)


def comb_stack(inp, period, plots, template, start_offset, clip):
    """
    Remove periodic glitches from signal

    Input:
        inp (obspy.trace)  : the series to be deglitched
        period (float) : the glitch period (in seconds)
        plots (boolean) : produce plots
        dt (float)      : sample interval (seconds)
        template (obspy.trace) : ones and zeros indicating which parts of inp
                                 are to be used
        start_offset (float) : start the dirac comb at this offset from the
                           start_offset (places the glitch well within each
                           "slice")
        clip (tuple):   low and high values to clip trace at when calculating
                      average glitch

    Output:
        out (obspy.trace) : cleaned series (inp - y)
        y   (obspy.trace) : the glitch which was subtracted from inp
                            (conv(comb,x))
        x   (np.array)   :  averaged glitch
        xm  (np.array)   :  x starting one sample earlier
        xp  (np.array)   :  x starting one sample later
        c  (np.array)    :  comb
        ng  (float)      :  number of glitches used to make average
    """
    print(f'\tRunning comb_stack({period=}s)', flush=True)
    dt = inp.stats.delta
    samps_per_period = period/dt
    n = inp.stats.npts
    rp = int(np.round(samps_per_period))

    # -----------------------------------------
    # Create the Dirac comb, starting at start_offset
    if DEBUG:
        print('COMB_STACK(): Creating Dirac comb')
    if start_offset == 0:
        n_pulses, off, c = comb(n, samps_per_period)
        offset_samps = 0
    else:
        offset_samps = M.floor(start_offset/dt)
        n_pulses, off, c = comb(n-offset_samps, samps_per_period)
        c = np.hstack((np.zeros(offset_samps), c))

    # Create the broken Dirac comb (teeth missing for slices containing
    # template==0
    if DEBUG:
        print('COMB_STACK(): Creating broken Dirac comb...',
              flush=True, end='')
        tic = time.time()
    cb = c.copy()
    b_pulses = n_pulses
    for i in range(n_pulses):
        n1 = start_offset + M.floor(i*samps_per_period)
        n2 = n1 + rp
        if np.any(template.data[n1: n2] == 0):
            cb[n1: n2] = 0
            b_pulses = b_pulses-1
    if DEBUG:
        print('Took {:.1f} s'.format(time.time()-tic))

    if plots:
        plt.figure(102)
        plt.subplot(2, 1, 1)
        plt.plot(np.arange(c.size)*dt, c)
        plt.title('comb')
        plt.subplot(2, 1, 2)
        plt.plot(np.arange(cb.size)*dt, cb)
        plt.title('broken comb')
        plt.xlabel('Time (seconds)')
        plt.show()

    # -----------------------------------------
    # Calculate the average glitch by correlating
    # rp samples on either side of zero offset with the broken comb
    # ** Could make this faster by reducing the buffer on either
    # ** side to that which is extracted in the next section

    # adding x zeros to either side of the comb gives a length 2*rp+1 result
    # from correlate('valid')
    # buff=np.hstack( (np.zeros(rp),cb,np.zeros(rp)) )
    # Instead of buffering by rp (above), why not just buffer by a little more
    # than off (all we recover is ~-off-2:rp-off+2)
    nbuff = off + 5
    buff = np.hstack((np.zeros(rp), cb, np.zeros(nbuff)))
    if DEBUG:
        print('COMB_STACK(): Calculating average glitch by correlation...',
              flush=True, end='')
        tic = time.time()
    xx = np.correlate(inp.data. clip(clip[0], clip[1]), buff,
                      mode='valid') / b_pulses
    if DEBUG:
        print('Took {:.1f} s'.format(time.time()-tic))

    # Extract average glitch, as well as same
    # shifted one sample to right and one sample to left
    if DEBUG:
        print(f'{xx.size=}, {rp=}, {off=}')
        print('COMB_STACK(): Extracting average glitch')
#     x=xx[rp-off:2*rp-off+1]
#     xm=xx[rp-off-1:2*rp-off]
#     xp=xx[rp-off+1:2*rp-off+2]
#     xmm=xx[rp-off-2:2*rp-off-1]
#     xpp=xx[rp-off+2:2*rp-off+3]
    # Using nbuff instead of rp
    x = xx[nbuff-off: nbuff+rp-off+1]
    xm = xx[nbuff-off-1: nbuff+rp-off]
    xp = xx[nbuff-off+1: nbuff+rp-off+2]
    xmm = xx[nbuff-off-2: nbuff+rp-off-1]
    xpp = xx[nbuff-off+2: nbuff+rp-off+3]
    if plots:
        index = np.arange(rp) / (3600/dt)
        plt.figure(103)
        plt.plot(index, inp.data[offset_samps: offset_samps+rp], 'b',
                 label='signal')
        plt.plot(index, x[:rp], 'r', label='comb_stack', linewidth=2)
        plt.xlabel('hours')
        plt.title('blue: signal, red: comb_stack')
        plt.legend()
        plt.show()

    # Create synthetic glitch time series (using full comb)
    if DEBUG:
        print('COMB_STACK(): Creating synthetic glitch time series')
    y = inp.copy()
    yy = convolve(c, x)
    y.data = yy[off: off+n]

    if plots:
        index = np.arange(n) / (3600/dt)
        plt.figure(104)
        plt.plot(index, inp.data, 'b', label='signal')
        plt.plot(index, y.data, 'r', label='synglitch')
        plt.legend()
        plt.xlabel('hours')
        plt.title('Signal and synthetic glitch')
        plt.show()

    # Create corrected time series
    if DEBUG:
        print('COMB_STACK(): Creating corrected time series')
    out = inp.copy()
    out.data = inp.data - y.data
    return (out, y, x, xm, xp, c, b_pulses, xmm, xpp)


def comb(n, per):
    """
    Generate a Dirac comb with n/per impulses

    Set up in freq domain, then calculate the ifft

    WOULD PROBABLY BE MUCH FASTER IF CALCULATED FOR pow2(N), then cut down the
    RESULT TO LENGTH N

    :param n: length (samples)
    :param per: period (in samples)
    Output:
        n_pulses - number of pulses created
        off    - offset before pulses
        out    - Dirac comb
    """
    # print('Running comb(n={:g},per={:g})'.format(n,per),flush=True)
    if DEBUG:
        print('COMB(): Generating DIRAC comb')

    # CHANGE N TO A POWER OF TWO
    n_return = n
    n = 2**(M.ceil(M.log2(n_return)))
    if DEBUG:
        print('COMB(): n={:d}, npow2={:d}'.format(n_return, n))

    n_pulses = M.ceil(n/per)
    cl = np.zeros(n) + 1j*np.zeros(n)
    # pulses are shifted "off" samples to the right (to eliminate edge
    # effects?) to be removed after the convolution in comb_stack
    off = 48
    offex = 2*M.pi*1j*off/n
    fn2 = M.floor(n/2)
    fn4 = M.floor(n/4)
    if DEBUG:
        print('COMB(): Filling positive freq terms '
              '(0:{:d}) with fft of deltas...'.format(fn2), flush=True, end='')
        tic = time.time()

    # ORIGINAL METHOD
#     for l in range(fn2+1) :  # for l=0:fn2
#         lpn=l*per/n;
#         if np.abs(np.round(lpn)-lpn)<1e-9 :
#             cl[l]=n_pulses;
#         else :
#             q=np.exp(-2*M.pi*1j*lpn)
#             cl[l]=(q**n_pulses-1)/(q-1)
#         cl[l]=cl[l]*np.exp(-offex*l)
    # ALTERNATIVE (FASTER) METHOD
    lf = np.arange(fn2+1)
    lpn = lf*per/n
    q = np.exp(-2*M.pi*1j*lpn)
    iPulses = np.abs(np.round(lpn) - lpn) < 1e-9
    q[iPulses] = -1
    cl[:fn2+1] = np.exp(-offex*lf)*(q**n_pulses-1)/(q-1)
    cl[:fn2+1][iPulses] = n_pulses
    if DEBUG:
        print('Took {:.1f} seconds'.format(time.time() - tic))
        print('COMB(): multiplying upper half of frequencies '
              '({:d}:{:d}) by cos(0->pi)...'.format(fn4, fn2),
              flush=True, end='')
        tic = time.time()
    for ff in range(fn4, fn2+1):  # l=fn4:fn2
        cl[ff] *= (1 + M.cos((ff-fn4)/(fn2-fn4)*M.pi)) / 2.
    if fn2 == n/2:  # even # of elements
        # cl(fn2+2:n)=conj(cl(fn2:-1:2));
        cl[fn2+1: n] = np.conj(cl[fn2-1: 0: -1])
    else:
        # cl(fn2+2:n)=conj(cl(fn2+1:-1:2));
        cl[fn2+1:n] = np.conj(cl[fn2: 0: -1])
    if DEBUG:
        print('Took {:.1f} seconds'.format(time.time() - tic))
        print('COMB(): Calculating IFFT of {:d}-point series...'
              .format(cl.size), flush=True, end='')
        tic = time.time()

    out = np.real(ifft(cl))
    if DEBUG:
        print('Took {:.1f}s'.format(time.time() - tic))

    # CONVERT BACK TO ORIGINAL LENGTH
    n_pulses = M.ceil(n_return/per)
    out = out[: n_return]
    return n_pulses, off, out


def _remove_noisy(inp, template, start_offset, period, plot=False):
    """
    Sets template to zero for sections of the data with anomalously high noise

    Affects indices corresponding to the middle half of a slice (based on
    start_offset and period), to allow for subsequent small changes in period

    Inputs:
        inp (obspy.trace): the original data
        template (obspy.trace): the existing template (may already have zeros)
        start_offset (float): the offset (seconds) of the first slice
        period (float): the interval at which the data will be cut into slices
    Output:
        template (obspy.trace): the new template
  """
    dt = inp.stats.delta
    samp_per = period / dt
    rp = np.round(samp_per)
    n = inp.stats.npts - M.floor(start_offset/dt)
    n_slices = M.ceil(n / samp_per)

    # Calculate the variance of each data slice
    slices = list()
    for i in range(n_slices):
        n1 = start_offset + M.floor(i*samp_per)
        n2 = n1+rp
        # Don't go past end
        if n2 > n+1:
            continue
        # Ignore the first & last 10 samples to accomodate small changes in
        # period
        var = sum((inp.data[n1+10: n2-10] * template.data[n1+10: n2-10])**2)
        hasZeros = np.any(template.data[n1+10:n2-10] == 0)
        slices.append({'var': var, "n1": n1, "n2": n2, 'hasZeros': hasZeros})

    # Calculate the median and std of variances of slices without zeros
    median = np.median([x['var'] for x in slices if x['hasZeros'] is False])
    sigma = np.std([x['var'] for x in slices if x['hasZeros'] is False])

    # For slices with abnormally high variance, set template indices to zero
    print('\tRejecting high-variance slices (>{:.1e}+3*{:.1e})...'
          .format(var, sigma), end='')
    nRejected = 0
    for slice in slices:
        if slice['var'] <= median + 3*sigma:
            slice['muted'] = False
        else:
            nRejected += 1
            slice['muted'] = True
            n1 = slice['n1']
            n2 = slice['n2']
            quarter_slice = M.floor((n2-n1)/4)
            nStart = n1 + quarter_slice
            nEnd = n1 + 3*quarter_slice
            template.data[nStart:nEnd] = 0
    if nRejected:
        print('{:d} of {:d} rejected'.format(nRejected, len(slices)))
    else:
        print('none rejected')

    # Plot results
    if plot:
        stack = stack_data(inp.data[M.floor(start_offset/dt):], period)
        plt.figure(200)
        nrows, ncols = stack.shape
        time = np.arange(nrows)*dt
        hasZeros = [i for i in range(ncols) if slices[i]['hasZeros'] is True]
        isMuted = [i for i in range(ncols) if slices[i]['muted'] is True]
        isUsed = [i for i in range(ncols)
                  if slices[i]['hasZeroes'] is False
                  and slices[i]['muted'] is False]
        plt.plot(time, stack[:, hasZeros], 'r:', label='has Zeros')
        plt.plot(time, stack[:, isMuted], 'g--', label='high Variance (muted)')
        plt.plot(time, stack[:, isUsed], 'b-',
                 label='used for glitch calculation')

    return template
