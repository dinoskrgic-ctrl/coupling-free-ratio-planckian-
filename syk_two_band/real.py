"""Real-frequency solver for the same model (r = 0: frozen SYK bath), and the Kubo transport rate.

Spectral self-energies (rho_S = -Im S^R / pi):
  rho_Sf(w) = Uf^2 int A_f A_f A_f delta(w-w1-w2+w3) [(1-f1)(1-f2)f3 + f1 f2 (1-f3)]
  rho_Sc(w) = U^2  int A_c(w1) A_f(w2) A_f(w3) delta(...) [same occupation factors]
evaluated as products in time on a periodic grid; Re S^R by the discrete Hilbert transform (theta(t) in time).
"""
import numpy as np
from matsu import g_semi

def fermi(x, T):
    return 0.5 * (1 - np.tanh(x / (2 * T)))

class RGrid:
    def __init__(self, L=4.0, N=2**21):
        self.N = N; self.L = L
        self.dw = 2 * L / N
        self.w = -L + self.dw * np.arange(N)
        self.dt = 2 * np.pi / (N * self.dw)
        k = np.arange(N)
        self.ph = np.exp(-1j * self.w[0] * self.dt * k)      # e^{-i w_min t_k}
        th = np.zeros(N); th[1:N // 2] = 1.0; th[0] = 0.5
        self.theta = th

    def to_t(self, X):
        """x(t_k) = int dw X(w) e^{-i w t_k}"""
        return self.dw * self.ph * np.fft.fft(X)

    def to_w(self, x):
        """X(w_j) = (1/2pi) int dt e^{i w t} x(t)  (inverse of to_t)"""
        return (self.dt / (2 * np.pi)) * self.N * np.fft.ifft(x * np.conj(self.ph))

    def retarded(self, rho):
        """S^R(w) = int dx rho(x)/(w - x + i0) via S^R(t) = -i theta(t) int rho e^{-ixt}"""
        r = self.to_t(rho)
        St = -1j * self.theta * r
        return 2 * np.pi * self.to_w(St)

def rho_sigma(g, Aa, Ab, T, coupling):
    """coupling^2 * int Aa(w1) Ab(w2) Ab(w3) d(w-w1-w2+w3)[(1-f1)(1-f2)f3 + f1 f2 (1-f3)]"""
    f = fermi(g.w, T)
    a_g = g.to_t(Aa * (1 - f)); a_l = g.to_t(Aa * f)
    b_g = g.to_t(Ab * (1 - f)); b_l = g.to_t(Ab * f)
    P = a_g * b_g * np.conj(b_l) + a_l * b_l * np.conj(b_g)
    rho = coupling**2 * g.to_w(P).real
    return np.clip(rho, 0, None)

def solve_f(g, T, Uf, A0=None, tol=1e-9, maxit=3000, mix=0.3, verbose=False):
    w = g.w
    if A0 is None:
        A = (1 / np.pi) * (Uf / (w**2 + Uf**2))
    else:
        A = A0
    for it in range(maxit):
        rho = rho_sigma(g, A, A, T, Uf)
        S = g.retarded(rho)
        S = 1j * S.imag * 0 + S       # keep full
        G = 1 / (w - S + 1e-14j)
        An = -G.imag / np.pi
        err = np.max(np.abs(An - A)) / np.max(A)
        A = (1 - mix) * A + mix * An
        if verbose and it % 50 == 0: print('f', it, err)
        if err < tol: break
    return A, S, it, err

def solve_c(g, T, D, U, Af, tol=1e-10, maxit=500, mix=0.7):
    w = g.w
    Gc = g_semi(w + 1e-12j, D)
    Ac = -Gc.imag / np.pi
    for it in range(maxit):
        rho = rho_sigma(g, Ac, Af, T, U)
        S = g.retarded(rho)
        Gn = g_semi(w - S + 1e-12j, D)
        An = -Gn.imag / np.pi
        err = np.max(np.abs(An - Ac)) / np.max(Ac)
        Ac = (1 - mix) * Ac + mix * An
        if err < tol: break
    return Ac, S, it, err

def tau_tr(g, T, D, S):
    """Kubo transport time without vertex corrections (random couplings):
    tau_tr = pi int dw (-f') int de Phi(e)/Phi(0) A(e,w)^2, Phi ~ (D^2-e^2)^{3/2} (Bethe-lattice transport function).
    For Gamma = -Im S << D the e integral is Phi(x)/(2 pi Gamma) with x = w - Re S, up to O(Gamma^2/D^2)."""
    w = g.w; m = np.abs(w) < 40 * T
    x = w[m] - S[m].real; Gam = -S[m].imag
    Phi = np.clip(1 - (x / D)**2, 0, None) ** 1.5
    mf = 1 / (4 * T * np.cosh(w[m] / (2 * T))**2)
    return np.trapezoid(mf * Phi / (2 * Gam), w[m])

def solve_fc(g, T, D, U, Uf, r, Af0=None, tol=1e-9, maxit=3000, mix=0.3):
    """Joint solution with bath feedback: Sigma_f = SYK + r U^2 [A_f A_c A_c] term (r = N_c/N_f)."""
    w = g.w
    Af = Af0 if Af0 is not None else (1 / np.pi) * (Uf / (w**2 + Uf**2))
    Ac = -g_semi(w + 1e-12j, D).imag / np.pi
    for it in range(maxit):
        rf = rho_sigma(g, Af, Af, T, Uf) + rho_sigma(g, Af, Ac, T, np.sqrt(r) * U)
        Sf = g.retarded(rf)
        Afn = -(1 / (w - Sf + 1e-14j)).imag / np.pi
        rc = rho_sigma(g, Ac, Af, T, U)
        Sc = g.retarded(rc)
        Acn = -g_semi(w - Sc + 1e-12j, D).imag / np.pi
        err = max(np.max(np.abs(Afn - Af)) / np.max(Af), np.max(np.abs(Acn - Ac)) / np.max(Ac))
        Af = (1 - mix) * Af + mix * Afn; Ac = (1 - mix) * Ac + mix * Acn
        if err < tol: break
    return Af, Ac, Sf, Sc, it, err
