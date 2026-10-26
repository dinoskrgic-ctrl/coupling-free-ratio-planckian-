"""Bath functionals shared by all v6 scripts (Supplemental Material Secs. S2-S3).

A local marginal bath has Eliashberg function a2F(W,T) = (lam/2) w(W/T) with w -> 1 for W >> T.
  S0(w)  = beta = 2 int_0^inf w(x)/sinh(x) dx           : -Im Sigma(0,T) = (pi/2) lam beta T
  Cb(w)  = transport coefficient, alpha = Cb*lam         : hbar/tau_tr = alpha kB T
  optical_model : extended-Drude hbar/tau(omega)/kT with the full Kramers-Kronig self-energy
Units hbar = kB = 1, T = 1 unless stated."""
import numpy as np
from functools import lru_cache
from scipy.integrate import quad
from scipy.special import zeta
from scipy.signal import hilbert

QS = 4*np.pi**3/(21*zeta(3))            # conformal Q = 4.913
CB_CONF = 2*np.pi**3/(7*zeta(3))        # conformal C_b = 7.370
w_conf = lambda x: np.tanh(x/2)


def w_arc(r):
    return lambda x: (2/np.pi)*np.arctan(x/r)


def S0(w):
    return 2*quad(lambda x: w(x)/np.sinh(x), 1e-9, 60, limit=400)[0]


def _rate(e, w):
    f = lambda z: 0.5*(1-np.tanh(z/2))
    return quad(lambda W: w(W)*(2/np.expm1(W) + f(W+e) + f(W-e)), 1e-9, abs(e)+60, limit=400)[0]


def _dc(w):
    xs = np.linspace(-25, 25, 201); wt = 1/(4*np.cosh(xs/2)**2)
    return 1/np.trapezoid(wt/np.array([_rate(x, w) for x in xs]), xs)


@lru_cache(maxsize=1)
def _dc_conf():
    return _dc(w_conf)


def Cb(w):
    return CB_CONF*_dc(w)/_dc_conf()


@lru_cache(maxsize=1)
def arctan_table():
    """(r, beta, C_b) on a grid for w = (2/pi) arctan(x/r); beta and C_b fall monotonically with r."""
    rs = np.logspace(np.log10(0.15), np.log10(8.0), 48)
    be = np.array([S0(w_arc(r)) for r in rs]); cb = np.array([Cb(w_arc(r)) for r in rs])
    return rs, be, cb


def beta_of_Cb(c):
    rs, be, cb = arctan_table(); return np.interp(c, cb[::-1], be[::-1])


def r_of_Cb(c):
    rs, be, cb = arctan_table(); return np.exp(np.interp(c, cb[::-1], np.log(rs[::-1])))


def Cb_of_beta(b):
    rs, be, cb = arctan_table(); return np.interp(b, be[::-1], cb[::-1])


def optical_model(w, lam, xs, wc=400.0, emax=700.0, de=0.1):
    """Extended-Drude hbar/tau(omega)/kT for a2F = (lam/2) w(W/T), T = 1, flat up to the cutoff wc (in units of T)."""
    e = np.arange(-emax, emax+de/2, de)
    W = np.linspace(1e-6, wc, 6000); dW = W[1]-W[0]; a2f = 0.5*lam*w(W); nb = 1/np.expm1(W)
    f = lambda x: 0.5*(1-np.tanh(x/2))
    ims = np.array([-np.pi*np.sum(a2f*(2*nb+f(W+x)+f(W-x)))*dW for x in e])
    res = np.imag(hilbert(ims-ims[0]))       # Kramers-Kronig on the symmetric grid
    S = res + 1j*ims
    fe = f(e); out = []
    for x in xs:
        k = max(1, int(round(x/de))); w_ = k*de
        sig = (1j/w_)*np.sum((fe[:-k]-fe[k:])/(w_+np.conj(S[:-k])-S[k:]))*de
        out.append((1j/sig - w_).imag)
    return np.array(out)
