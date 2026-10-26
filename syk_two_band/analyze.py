"""Extract B from the free energy and alpha from the Kubo rate; form Q = alpha gamma_b / B."""
import numpy as np, glob, sys
from scipy.special import zeta
base = ''+__import__('os').path.dirname(__import__('os').path.abspath(__file__))+f'/out/'
Qs = 4 * np.pi**3 / (21 * zeta(3)); Cb = 2 * np.pi**3 / (7 * zeta(3))
gb = 2 * np.pi / 3          # (pi^2/3) rho0(0), semicircle D = 1, per flavor

def fitB(T, F, Tmax, extra=True):
    m = T <= Tmax
    T, F = T[m], F[m]
    cols = [np.ones_like(T), T**2, T**2 * np.log(T)]
    if extra: cols += [T**3, T**3 * np.log(T)]
    A = np.array(cols).T
    sc = np.max(np.abs(A), axis=0); A2 = A / sc
    c, *_ = np.linalg.lstsq(A2, F, rcond=None); c = c / sc
    res = F - A @ c
    return 2 * c[2], np.sqrt(np.mean(res**2)), c

if __name__ == '__main__':
    for f in sorted(glob.glob(base + 'F_U*_dt*.npy')):
        d = np.load(f); T, F = d[:, 0], d[:, 1]
        for Tmax in (1e-3, 2e-3, 4e-3):
            for extra in (False, True):
                B, rms, c = fitB(T, F, Tmax, extra)
                print(f"{f.split('/')[-1]:22s} Tmax={Tmax:.0e} T3terms={extra!s:5s}  B={B:.5f}  lam_C=(2/3)B/gb={2*B/(3*gb):.5f}  rms={rms:.1e}")
