"""
Python version of matlab compliance functions available on Wayne Crawford's
website as of October 2017

Authors:  A. Doran, W Crawford
"""
import math
import numpy as np


def gravd(W,h):
    """
    Linear ocean surface gravity wave dispersion
    
    Args:
        W (:class:`numpy.ndarray`): angular frequencies (rad/s)
        h (float): water depth (m)

    Returns:
        K (:class:`numpy.ndarray`): wavenumbers (rad/m)
    """
    # W must be array
    if not insintance(W, np.ndarray):
        W = np.array([W])
    G = 9.79329
    N = len(W)
    W2 = W*W
    kDEEP = W2/G
    kSHAL = W/(np.sqrt(G*h))
    erDEEP = np.ones(np.shape(W)) - G*kDEEP*_dtanh(kDEEP*h)/W2
    one = np.ones(np.shape(W))
    d = np.copy(one)
    done = np.zeros(np.shape(W))
    nd = np.where(done==0)
    
    k1 = np.copy(kDEEP)
    k2 = np.copy(kSHAL)
    e1 = np.copy(erDEEP)
    ktemp = np.copy(done)
    e2 = np.copy(done)
        
    while True:
        e2[nd] = one[nd] - G*k2[nd] * _dtanh(k2[nd]*h)/W2[nd]
        d = e2*e2
        done = d<1e-20
        if done.all():
            K = k2
            break
        
        nd=np.where(done==0)
        ktemp[nd] = k1[nd]-e1[nd]*(k2[nd]-k1[nd])/(e2[nd]-e1[nd])
        k1[nd] = k2[nd]
        k2[nd] = ktemp[nd]
        e1[nd] = e2[nd]
    
    return K


def raydep(P,om,d,ro,vp2,vs2):
    """
    Propagator matrix solutionn for P-SV waves, minor vector method
    
    Args: 
        P (float): surface wave slowness (s/m)
        om (float): surface wave angular frequency (radians/sec)
        d (:class:`numpy.ndarray`): thicknesses of the model layers (meters?)
        rho (:class:`numpy.ndarray`): density of the layer (kg/m^3) (= gm/cc * 1000)
        vp2 (:class:`numpy.ndarray`): compressional velocity squared (m/s)^2
        vs2 (:class:`numpy.ndarray`): shear velocity squared (m/s)^2
        
    
    Returns:
        (list):
            v (:class:`numpy.ndarray`): vertical velocity AT TOP OF EACH LAYER
            u (:class:`numpy.ndarray`): horizontal velocity AT TOP OF EACH LAYER
            zz (:class:`numpy.ndarray`): vertical stress AT TOP OF EACH LAYER
            zx (:class:`numpy.ndarray`): horizontal stress AT TOP OF EACH LAYER
   
    Notes:
        d, rho, vp2 and vs2 must be same length
        (Normalized compliance = -k*v/(omega*sigzz) )
    """
    mu = ro * vs2;
    n = len(d);
    ist = n-1;
    ysav = 0;
    psq = P*P;
    r2 = 2*mu[ist]*P;
    #% R and S are the "Wavenumbers" of compress and shear waves in botlayer
    #% RoW and SoW are divided by ang freq
    RoW = np.sqrt(psq- 1/vp2[ist]);
    SoW = np.sqrt(psq- 1/vs2[ist]);
    ym = np.zeros((ist+1,5));
    i = ist;
    y=np.zeros((5,))    # Minor vector matrix
    x=np.zeros((i+1,4)) # Stress-displacement vector:
                        #      (vert vel, horiz vel, vert stress, hor stress)
    y[2] =  RoW
    y[3] = -SoW
    y[0] = (RoW*SoW - psq) / ro[i]
    y[1] = r2*y[0] + P
    y[4] = ro[i] - r2*(P + y[1])
    ym[i,:] = y
    # *****PROPAGATE UP LAYERS*********
    while i > 0:
        i = i-1
        ha = psq - 1/vp2[i]
        ca,sa = _argdtray(om*d[i], ha)
        hb = psq - 1/vs2[i]
        cb,sb = _argdtray(om*d[i], hb)
        hbs = hb*sb
        has = ha*sa
        r1 = 1 /r o[i]
        r2 = 2 * mu[i] * P
        b1 = r2*y[0] - y[1]
        g3 = ( y[4] + r2*(y[1]-b1) ) * r1
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
        ym[i,:] = y
    
    de = y[5-1]/np.sqrt(y[1-1]*y[1-1] + y[2-1]*y[2-1])
    ynorm = 1/y[3-1]
    y[1-1:4] = np.array([0 , -ynorm,  0,  0])
    # *****PROPAGATE BACK DOWN LAYERS*********
    while i <= ist:
        x[i, 0] = -ym[i, 1]*y[1-1] - ym[i,3-1]*y[2-1] + ym[i,1-1]*y[4-1]
        x[i, 1] = -ym[i, 3]*y[1-1] + ym[i,2-1]*y[2-1] - ym[i,1-1]*y[3-1]
        x[i, 2] = -ym[i, 4]*y[2-1] - ym[i,2-1]*y[3-1] - ym[i,4-1]*y[4-1]
        x[i, 3] =  ym[i, 4]*y[1-1] - ym[i,3-1]*y[3-1] + ym[i,2-1]*y[4-1]
        ls = i;
        if i >=2-1:
            sum = abs( x[i,1-1] + i*x[i,2-1]);
            pbsq = 1/vs2[i];
            if sum < 1e-4:
                break
                
        ha = psq - 1/vp2[i];
        ca, sa = _argdtray(om*d[i],ha);
        hb = psq-1/vs2[i];
        cb,sb = _argdtray(om*d[i],hb);
        hbs = hb*sb;
        has = ha*sa;
        r2 = 2*P*mu[i];
        e2 = r2*y[2-1] - y[3-1];
        e3 = ro[i]*y[2-1] - P*e2;
        e4 = r2*y[1-1] - y[4-1];
        e1 = ro[i]*y[1-1] - P*e4;
        e6 = ca*e2 - sa*e1;
        e8 = cb*e4 - sb*e3;
        y[1-1] = (ca*e1 - has*e2+P*e8) / ro[i];
        y[2-1] = (cb*e3 - hbs*e4+P*e6) / ro[i];
        y[3-1] = r2*y[2-1] - e6;
        y[4-1] = r2*y[1-1] - e8;
        i = i+1;
    #
    #if x(1,3) == 0
    #  error('vertical surface stress = 0 in DETRAY');
    #end
    ist = ls;
    v =	 x[:,1-1];
    u =  x[:,2-1];
    zz = x[:,3-1];
    zx = x[:,4-1];
    
    return v, u, zz, zx



def calc_norm_compliance(depth,freq,model):
    """ calculate normalized compliance for a model and water depth
    
    model=[thick(m)  rho(g/cc) vp(m/s) vs(m/s)]
    """

    thick=model[:,0];
    rho=model[:,1];
    vpsq=model[:,2]*model[:,2];
    vssq=model[:,3]*model[:,3];
    omega = 2*np.pi*freq;
    k = gravd(omega,depth);
    p = k / omega;
    
    ncomp=np.zeros((len(p)))
    
    for i in np.arange((len(p))):
        v, u, sigzz, sigzx = raydep( p[i],omega[i],thick,rho,vpsq,vssq)
        ncomp[i] = -k[i]*v[1-1]/(omega[i]*sigzz[1-1]);
    
    return ncomp


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


def _argdtray(wd,h):
    hh = np.sqrt(abs(h))    # magnitude of wavenumber/freq
    th = wd * hh 			# number of waves (or e-foldings) in layer in radians
    if th >= 1.5e-14:
        if h <= 0:          #  propagating wave
            c =  np.cos(th)
            s = -np.sin(th) / hh
        else:			    # evenescent wave
           	d=np.exp(th);
           	c =  0.5*(d + 1/d)
           	s = -0.5*(d - 1/d)/hh
    else:
        c = 1
        s = -wd
    
    return c,s



