"""Imaginary-time large-N solver for the two-band SYK model at half filling.

c fermions: semicircular band of half-width D (rho0 = 2 sqrt(D^2-e^2)/(pi D^2)), N_c flavors.
f fermions: SYK_4 islands, coupling Uf, N_f flavors.
Inter-band:  H_cf = sum U_{ij;kl} c_i^dag c_j f_k^dag f_l, random, variance U^2 (per c flavor, summed over f).

Large N:  Sigma_c(tau) = -U^2 G_c(tau) G_f(tau) G_f(-tau)
          Sigma_f(tau) = -Uf^2 G_f(tau)^2 G_f(-tau)  - r U^2 G_f(tau) G_c(tau) G_c(-tau),  r = N_c/N_f
With particle-hole symmetry G(-tau) = -G(tau) for 0<tau<beta.
r = 0 is the frozen-bath limit (N_f >> N_c); r = 1 is the original equal-flavor model.

Free energy per c flavor, O(N_c) part, relative to decoupled bands (Luttinger-Ward):
  dF = -T sum_n [L(iw-Sc) - L(iw)] - T sum_n Sc Gc + Phi_cf
       + r * ( f-sector change: -T sum_n [ln(iw-Sf) - ln(iw-Sf0)] - T sum (Sf Gf - Sf0 Gf0) + Phi_ff - Phi_ff0 )
  with Phi_cf = -(U^2/2) int_0^beta Gc(t)Gc(-t)Gf(t)Gf(-t) dt,   Phi_ff = -(Uf^2/4) int Gf(t)^2 Gf(-t)^2 dt.
For r = 0 the f terms vanish identically (f frozen, its own free energy is O(N_f) and coupling-independent).
"""
import numpy as np

def L_semi(z, D):
    """int rho0(e) ln(z-e) de for the semicircle; branch with s ~ z at infinity."""
    s = np.sqrt(z - D) * np.sqrt(z + D)
    return (z * z - z * s) / D**2 + np.log((z + s) / 2) - 0.5

def g_semi(z, D):
    s = np.sqrt(z - D) * np.sqrt(z + D)
    return 2 * (z - s) / D**2

class Grid:
    def __init__(self, T, dtau):
        self.T = T; self.beta = 1.0 / T
        N = int(round(self.beta / dtau)); N += N % 2
        self.N = N; self.dt = self.beta / N
        n = np.fft.fftfreq(N, d=1.0 / N).astype(np.int64)     # n index in fft order
        self.n = n
        self.w = (2 * n + 1) * np.pi * T                         # fermionic Matsubara, fft order
        self.tau = np.arange(N) * self.dt
        th = self.w * self.dt
        # Filon (piecewise-linear) weights: phi = int_0^1 (1-u) e^{i th u}, psi = int_0^1 u e^{i th u}
        small = np.abs(th) < 1e-4
        ths = np.where(small, 1.0, th)
        e = np.exp(1j * ths)
        psi = e / (1j * ths) - (e - 1) / (1j * ths) ** 2
        phi = (e - 1) / (1j * ths) - psi
        phi = np.where(small, 0.5 + 1j * th / 6, phi)
        psi = np.where(small, 0.5 + 1j * th / 3, psi)
        self.phi, self.psi = phi, psi
        self.ph_t = np.exp(-1j * np.pi * np.arange(N) / N)

    def w2t(self, Gw, M2):
        """G(tau_k), 0<=tau<beta, from G(iw_n), subtracting 1/iw and M2/(iw)^3 tails (real, PH symmetric)."""
        iw = 1j * self.w
        x = Gw - 1 / iw - M2 / iw**3
        G = self.T * self.ph_t * np.fft.fft(x)
        tau = self.tau
        return G.real - 0.5 + M2 * tau * (self.beta - tau) / 4

    def t2w(self, St, St_beta):
        """int_0^beta S(tau) e^{i w tau} dtau, S piecewise linear on the grid; St_beta = S(beta^-)."""
        N = self.N
        s = N * np.fft.ifft(St * np.conj(self.ph_t))            # sum_j S_j e^{i w tau_j}
        emi = np.exp(-1j * self.w * self.dt)
        return self.dt * ((self.phi + self.psi * emi) * s - self.psi * emi * St[0] - St_beta * self.psi * emi)

def solve(T, D=1.0, Uf=0.05, U=0.5, r=0.0, dtau=0.02, tol=1e-12, maxit=4000, mix=0.5, init=None, verbose=False):
    gr = Grid(T, dtau)
    iw = 1j * gr.w
    # f: start from conformal-ish guess
    if init is None:
        Gf = 1 / (iw - (-1j * np.sign(gr.w) * np.sqrt(Uf * np.abs(gr.w)) * 0.8))
        Gc = g_semi(iw, D)
    else:
        Gf, Gc = init
    m1f = -2 * (Uf**2) * (-0.5)**3 * 0      # initial tail guesses (updated each iteration)
    m1c = 0.0
    inner = np.abs(gr.n) < gr.N // 4
    for it in range(maxit):
        gf = gr.w2t(Gf, m1f); gc = gr.w2t(Gc, D**2 / 4 + m1c)
        Sf_t = Uf**2 * gf**3 + r * U**2 * gf * gc**2
        Sc_t = U**2 * gc * gf**2
        m1f = -2 * Sf_t[0]; m1c = -2 * Sc_t[0]
        Sf = 1j * gr.t2w(Sf_t, Sf_t[0]).imag; Sc = 1j * gr.t2w(Sc_t, Sc_t[0]).imag
        Gf_new = 1 / (iw - Sf)
        Gc_new = g_semi(iw - Sc, D)
        err = max(np.max(np.abs(Gf_new - Gf)[inner]), np.max(np.abs(Gc_new - Gc)[inner]))
        Gf = (1 - mix) * Gf + mix * Gf_new
        Gc = (1 - mix) * Gc + mix * Gc_new
        if verbose and it % 50 == 0: print(it, err)
        if err < tol: break
    gf = gr.w2t(Gf, m1f); gc = gr.w2t(Gc, D**2 / 4 + m1c)
    Sf_t = Uf**2 * gf**3 + r * U**2 * gf * gc**2
    Sc_t = U**2 * gc * gf**2
    Sf = 1j * gr.t2w(Sf_t, Sf_t[0]).imag; Sc = 1j * gr.t2w(Sc_t, Sc_t[0]).imag
    return dict(gr=gr, Gf=Gf, Gc=Gc, gf=gf, gc=gc, Sf=Sf, Sc=Sc, it=it, err=err, D=D, Uf=Uf, U=U, r=r)

def trap(y, dt):
    """trapezoid over [0,beta] with y(beta) = y(0) (PH symmetric integrands)."""
    return dt * (np.sum(y) - 0.5 * y[0] + 0.5 * y[0])

def free_energy(sol, sol0_f=None):
    """O(N_c) free energy per c flavor relative to decoupled bands."""
    gr, D, U, Uf, r = sol['gr'], sol['D'], sol['U'], sol['Uf'], sol['r']
    T = gr.T; iw = 1j * gr.w
    Sc, Gc, gc, gf = sol['Sc'], sol['Gc'], sol['gc'], sol['gf']
    term1 = -T * np.sum(L_semi(iw - Sc, D) - L_semi(iw, D)).real
    term2 = -T * np.sum(Sc * Gc).real
    Phi = -(U**2 / 2) * trap(gc**2 * gf**2, gr.dt)
    dF = term1 + term2 + Phi
    if r != 0:
        Sf, Gf = sol['Sf'], sol['Gf']; Sf0, Gf0, gf0 = sol0_f
        t1 = -T * np.sum(np.log((iw - Sf) / (iw - Sf0))).real
        t2 = -T * np.sum(Sf * Gf - Sf0 * Gf0).real
        Pff = -(Uf**2 / 4) * (trap(gf**4, gr.dt) - trap(gf0**4, gr.dt))
        # Phi_cf already counts the c-f term once per c flavor; with r != 0 the same functional is shared
        dF += r * (t1 + t2 + Pff)
    return dF
