"""
Python compliance functions:
- gravd
- raydep
- zp_to_norm_compliance
- frf_to_compliance
- calc_norm_compliance

Authors:  A. Doran, W Crawford
"""
from copy import deepcopy

# import math
import numpy as np

# from .earth_model import EarthModel1D


def gravd(W, h):
    """
    Linear ocean surface gravity wave dispersion

    Args:
        W (:class:`numpy.ndarray`): angular frequencies (rad/s)
        h (float): water depth (m)

    Returns:
        K (:class:`numpy.ndarray`): wavenumbers (rad/m)

    >>> f = np.array([0.0001, 0.001, 0.01, 0.1, 1, 10])
    >>> K = gravd(2 * np.pi * f, 2000)
    >>> wlen = 2 * np.pi * np.power(K, -1)
    >>> np.set_printoptions(precision=1)
    >>> print(f'{K=} rad/m')
    K=array([4.5e-06, 4.5e-05, 5.2e-04, 4.0e-02, 4.0e+00, 4.0e+02]) rad/m
    >>> print(f'{wlen=} m')
    wlen=array([1.4e+06, 1.4e+05, 1.2e+04, 1.6e+02, 1.6e+00, 1.6e-02]) m
    >>> print(f'c={f*wlen} m/s')
    c=[140.  139.8 121.1  15.6   1.6   0.2] m/s
    """
    # W must be array
    if not isinstance(W, np.ndarray):
        W = np.array([W])
    if np.any(W < 0):
        raise ValueError('there are omegas <= 0')
    G = 9.79329
    # N = len(W)
    W2 = W*W
    kDEEP = W2/G
    kSHAL = W/(np.sqrt(G*h))
    erDEEP = np.ones(np.shape(W)) - G*kDEEP*_dtanh(kDEEP*h)/W2
    one = np.ones(np.shape(W))
    d = np.copy(one)
    done = np.zeros(np.shape(W))
    done[W == 0] = 1   # if W==0, k is also zero
    nd = np.where(done == 0)

    k1 = np.copy(kDEEP)
    k2 = np.copy(kSHAL)
    e1 = np.copy(erDEEP)
    ktemp = np.copy(done)
    e2 = np.copy(done)
    e2[W == 0] = 0

    while True:
        e2[nd] = one[nd] - G*k2[nd] * _dtanh(k2[nd]*h)/W2[nd]
        d = e2*e2
        done = d < 1e-20
        if done.all():
            K = k2
            break
        nd = np.where(done == 0)
        ktemp[nd] = k1[nd]-e1[nd]*(k2[nd]-k1[nd])/(e2[nd]-e1[nd])
        k1[nd] = k2[nd]
        k2[nd] = ktemp[nd]
        e1[nd] = e2[nd]
    return K


def raydep(P, om, d, ro, vp2, vs2):
    """
    Propagator matrix solutionn for P-SV waves, minor vector method

    Args:
        P (float): surface wave slowness (s/m)
        om (float): surface wave angular frequency (radians/sec)
        d (:class:`numpy.ndarray`): layer thicknesses (meters?)
        rho (:class:`numpy.ndarray`): layer densities (kg/m^3) (gm/cc * 1000)
        vp2 (:class:`numpy.ndarray`): layer P velocities squared (m/s)^2
        vs2 (:class:`numpy.ndarray`): layer shear velocities squared (m/s)^2

    Returns:
        (list): Parameters, each value is at layer top
            v (:class:`numpy.ndarray`): vertical velocity (m/s?)
            u (:class:`numpy.ndarray`): horizontal velocity (m/s?)
            zz (:class:`numpy.ndarray`): vertical stress (Pa?)
            zx (:class:`numpy.ndarray`): horizontal stress (Pa?)

    Notes:
        d, rho, vp2 and vs2 have one value for each layer (top to bottom),
            must be same length
        (Normalized compliance = -k*v/(omega*sigzz) )

    >>> P = 1/140    # Corresponds to 2000m depth, low freqs
    >>> om = 2 * np.pi * 0.005
    >>> d = np.array([1000, 1000, 1000, 3000, 3000])
    >>> rho = np.array([3000, 3000, 3000, 3000, 3000])
    >>> vp2 = np.array([3000**2, 4000**2, 5000**2, 7500**2, 8200**2])
    >>> vs2 = np.array([1600**2, 2300**2, 2800**2, 4300**2, 4700**2])
    >>> np.set_printoptions(precision=1)
    >>> raydep(P, om, d, rho, vp2, vs2)
    (array([1. , 0.7, 0.5, 0.4, 0.3]), array([ 1.8e-01,  6.0e-02,  1.2e-02,  5.2e-05, -5.8e-02]), array([-2.8e+08, -2.8e+08, -2.7e+08, -2.6e+08, -1.9e+08]), array([-0.0e+00,  3.1e+07,  5.4e+07,  7.6e+07,  1.0e+08]))
    """
    mu = ro * vs2
    n = len(d)
    ist = n-1
    # ysav = 0
    psq = P*P
    r2 = 2 * mu[ist] * P
    # R and S are the "Wavenumbers" of compress and shear waves in botlayer
    # RoW and SoW are divided by ang freq
    RoW = np.sqrt(psq - 1/vp2[ist])
    SoW = np.sqrt(psq - 1/vs2[ist])
    ym = np.zeros((ist+1, 5))
    i = ist
    y = np.zeros((5, ))     # Minor vector matrix
    # Stress-displacement vector: (vert vel, hor vel, vert stress, hor stress)
    x = np.zeros((i+1, 4))

    y[2] =  RoW
    y[3] = -SoW
    y[0] = (RoW*SoW - psq) / ro[i]
    y[1] = r2*y[0] + P
    y[4] = ro[i] - r2*(P + y[1])
    ym[i, :] = y
    # *****PROPAGATE UP LAYERS*********
    while i > 0:
        i = i-1
        ha = psq - 1/vp2[i]
        ca, sa = _argdtray(om*d[i], ha)
        hb = psq - 1/vs2[i]
        cb, sb = _argdtray(om*d[i], hb)
        hbs = hb*sb
        has = ha*sa
        r1 = 1 / ro[i]
        r2 = 2 * mu[i] * P
        b1 = r2*y[0] - y[1]
        g3 = (y[4] + r2*(y[1]-b1)) * r1
        g1 = b1 + P*g3
        g2 = ro[i]*y[0] - P*(g1+b1)
        e1 = cb*g2 - hbs*y[2]
        e2 = -sb*g2 + cb*y[2]
        e3 = cb*y[3] + hbs*g3
        e4 = sb*y[3] + cb*g3
        y[2] = ca*e2 - has*e4
        y[3] = sa*e1 + ca*e3
        g3 = ca*e4 - sa*e2
        b1 = g1 - P*g3
        y[0] = (ca*e1 + has*e3 + P*(g1+b1))*r1
        y[1] = r2*y[0] - b1
        y[4] = ro[i]*g3 - r2*(y[1] - b1)
        ym[i, :] = y

    # de = y[4]/np.sqrt(y[0]*y[0] + y[1]*y[1])
    ynorm = 1/y[2]
    y[0: 4] = np.array([0, -ynorm,  0,  0])
    # *****PROPAGATE BACK DOWN LAYERS*********
    while i <= ist:
        x[i, 0] = -ym[i, 1]*y[0] - ym[i, 2]*y[1] + ym[i, 0]*y[3]
        x[i, 1] = -ym[i, 3]*y[0] + ym[i, 1]*y[1] - ym[i, 0]*y[2]
        x[i, 2] = -ym[i, 4]*y[1] - ym[i, 1]*y[2] - ym[i, 3]*y[3]
        x[i, 3] =  ym[i, 4]*y[0] - ym[i, 2]*y[2] + ym[i, 1]*y[3]
        ls = i
        if i >= 1:
            sum = abs(x[i, 0] + i*x[i, 1])
            # pbsq = 1 / vs2[i]
            if sum < 1e-4:
                break

        ha = psq - 1/vp2[i]
        ca, sa = _argdtray(om*d[i], ha)
        hb = psq-1/vs2[i]
        cb, sb = _argdtray(om*d[i], hb)
        hbs = hb*sb
        has = ha*sa
        r2 = 2*P*mu[i]
        e2 = r2*y[1] - y[2]
        e3 = ro[i]*y[1] - P*e2
        e4 = r2*y[0] - y[3]
        e1 = ro[i]*y[0] - P*e4
        e6 = ca*e2 - sa*e1
        e8 = cb*e4 - sb*e3
        y[0] = (ca*e1 - has*e2+P*e8) / ro[i]
        y[1] = (cb*e3 - hbs*e4+P*e6) / ro[i]
        y[2] = r2*y[1] - e6
        y[3] = r2*y[0] - e8
        i = i+1
    #
    # if x(1,3) == 0
    #   error('vertical surface stress = 0 in DETRAY');
    # end
    ist = ls

    return x[:, 0], x[:, 1], x[:, 2], x[:, 3]


def compliance(depth, freq, model):
    """
    Calculate compliance of a model for a give water depth

    Args:
        depth (float): water depth (m)
        freq (:class:`numpy.nparray`): frequencies (1/s)
        model (:class:`EarthModel1D`): 1D earth model
    """
    if np.any(freq <= 0):
        raise ValueError('At least one freq <= 0: cannot calculate compliance')
    vpsq = model.vps * model.vps
    vssq = model.vss * model.vss
    omega = 2 * np.pi * freq
    k = gravd(omega, depth)
    ps = k / omega

    compl = np.zeros((len(ps)))
    for i in np.arange((len(ps))):
        v, _, sigzz, _ = raydep(ps[i], omega[i], model.thicks, model.rhos,
                                    vpsq, vssq)
        # If raydep returned complex values, would need to divide by a further
        # 1j to go from (m/s)/Pa to m/Pa.  Returned value should be
        # negative because seafloor is lowest (DOWN) under maxixum pressure,
        # for quasi-static
        compl[i] = v[0] / (omega[i] * sigzz[0])
    return compl


def calc_norm_compliance(depth, freq, model):
    """
    Calculate normalized compliance for a model and water depth
    
    norm compliance == k(omega) * compliance

    Args:
        depth (float): water depth (m)
        freq (:class:`numpy.nparray`): frequencies (1/s)
        model (:class:`EarthModel1D`): 1D earth model

    >>> depth = 2000
    >>> freqs = np.array([0.001, 0.003, 0.005, 0.01, 0.03])
    >>> model = EarthModel1D([[1000, 3000, 3000, 1600],
    ...                      [1000, 3000, 4000, 2300],
    ...                      [1000, 3000, 5000, 2800],
    ...                      [3000, 3000, 7500, 4300],
    ...                      [3000, 3000, 8200, 4700]])
    >>> np.set_printoptions(precision=1)
    >>> calc_norm_compliance(depth, freqs, model)
    array([-1.4e-11, -2.0e-11, -2.6e-11, -4.2e-11, -9.0e-11])
    """
    k = gravd(2 * np.pi * freq, depth)
    return k * compliance(depth, freq, model)

#     if np.any(freq <= 0):
#         raise ValueError('At least one freq <= 0: cannot calculate compliance')
#     vpsq = model.vps * model.vps
#     vssq = model.vss * model.vss
#     omega = 2 * np.pi * freq
#     k = gravd(omega, depth)
#     ps = k / omega
# 
#     ncomp = np.zeros((len(ps)))
#     for i in np.arange((len(ps))):
#         v, u, sigzz, sigzx = raydep(ps[i], omega[i], model.thicks, model.rhos,
#                                     vpsq, vssq)
#         ncomp[i] = (k[i] / omega[i]) * (v[0] / sigzz[0])
#     return ncomp


def zp_to_norm_compliance(freqs, zp, wdepth, z_units='M/S'):
    """
    Calculate normalized compliance from the z/p ratio, freqs and water depth

    normalized compliance is defined as k*Z/P, with Z in m and P in Pa.
    Its units are 1/Pa

    Args:
        freqs (:class:`numpy.nparray`): frequencies (1/s)
        zp (:class:`numpy.nparray`): vertical motion / pressure.  Pressure
            units are Pa, z_units are specified by the paramter z_units
        wdepth (float): water depth (m)
        z_units (str): z units, one of 'M', 'M/S' or 'M/S^2'
    """
    omega = 2 * np.pi * freqs
    k = gravd(omega, wdepth)
    if z_units.upper() == 'M':
        omega_term = np.ones(omega.shape())
    elif z_units.upper() == 'M/S':
        omega_term = omega**(-1)
    elif z_units.upper() == 'M/S^2':
        omega_term = omega**(-2)
    else:
        raise ValueError(f'Z_units ({z_units}) is not in ("M", "M/S", "M/S^2")')
    return zp * k * omega_term


def frf_to_compliance(xf, wdepth, z_units='M/S'):
    """
    Changes the response for each out_channel from z_units/Pa
    to 1/Pa (normalized compliance)

    Args:
        xf (:class:`tiskit.TransferFunction`): z/p transfer function(s)
        wdepth (float): water depth (m)
        z_units (str): z units, one of 'M', 'M/S' or 'M/S^2'
    """
    compl = deepcopy(xf)
    for oc in compl.output_channels:
        if not compl.output_units(oc).upper() == z_units:
            raise ValueError('output_units({}) ({}) != "{}"'.format(
                oc, compl.output_units(oc), z_units))
        if not compl.input_units.upper() == 'PA':
            raise ValueError(f'input_units ({compl.input_units}) != "PA"')
        orig_resp = compl.response(oc)
        new_resp = zp_to_norm_compliance(compl.freqs, orig_resp,
                                         wdepth, z_units)
        print(f'{new_resp/orig_resp=}')
        print(f'BEFORE {compl.response(oc)=}')
        compl.put_response(new_resp, oc)
        print(f'AFTER {compl.response(oc)=}')
        # compl._ds["response"].loc[dict(input=compl.input_channel,
        #                             output=oc)] = new_resp
    return compl


def _dtanh(x):
    """
    Stable hyperbolic tangent

    Args:
        x (:class:`numpy.ndarray`)
    """
    a = np.exp(x*(x <= 50))
    one = np.ones(np.shape(x))

    y = (abs(x) > 50) * (abs(x)/x) + (abs(x) <= 50)*((a-one/a) / (a+one/a))
    return y


def _argdtray(wd, h):
    hh = np.sqrt(abs(h))    # magnitude of wavenumber/freq
    th = wd * hh            # # of waves (or e-foldings) in layer (radians)
    if th >= 1.5e-14:
        if h <= 0:          # propagating wave
            c =  np.cos(th)
            s = -np.sin(th) / hh
        else:               # evenescent wave
            d=np.exp(th)
            c =  0.5*(d + 1/d)
            s = -0.5*(d - 1/d)/hh
    else:
        c = 1
        s = -wd
    return c, s


if __name__ == "__main__":
    import doctest
    doctest.testmod()
