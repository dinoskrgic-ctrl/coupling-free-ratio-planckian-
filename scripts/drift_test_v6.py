"""Drift test (Letter: 'The drift' and Fig. 2b; SM Sec. S8).

In the self-damped field r = m/(gamma T) falls logarithmically on cooling, so C_b(T), and with it the transport slope
alpha = C_b*lam, rises on cooling. With rho - rho_0 proportional to T*C_b(T), the local exponent is
n(T) = dln(rho - rho_0)/dlnT = 1 + dlnC_b/dlnT < 1. The conformal bath has n = 1 exactly.
This script (1) computes the predicted n(T) from the thermal-mass saddle (data/universal_curves.json) for bath
bandwidths E_L/gamma = 1 and 10 eV, and (2) measures n on three published resistivity curves:
  - Nd-LSCO p = 0.24, 33 T, 1-27 K           (Daou 2009, as replotted in Legros 2019 Fig. 1c)
  - LSCO p = 0.23, 48 T, 1.4-30 K            (Cooper 2009, as replotted in Legros 2019 Fig. 1b; 16 points)
  - LSCO p = 0.24, 0 T and 16 T, 20-264 K   (Michon 2023 Fig. 2a: the crystal used for the optics)
and (3) looks at the per-temperature zero-frequency optical rate (Michon 2023 Fig. 2b).
Writes ../outputs/drift_test.txt and ../data/drift_summary.json (used by make_figures_v6.py)."""
import os, json
import numpy as np
from scipy.optimize import curve_fit

H = os.path.dirname(os.path.abspath(__file__))
D = lambda f: os.path.join(H, '..', 'data', f)
KB = 8.617e-5                      # eV/K
lines = []
def say(s=''):
    print(s); lines.append(s)

# ---------------------------------------------------------------- (1) prediction
uc = json.load(open(D('universal_curves.json')))
CASES = {'Ec1.0_u0.3': 'u rho = 0.3', 'Ec1.0_u1': 'u rho = 1', 'Ec1.0_u1000': 'fixed amplitude (u -> inf)',
         'Ec0.3_u1000': 'fixed amplitude, E_c = 0.3 E_L', 'Ec3.0_u1000': 'fixed amplitude, E_c = 3 E_L'}
def n_pred(key, EL, TK):
    """local exponent 1 + dlnC_b/dlnT at temperatures TK (K) for bath bandwidth E_L/gamma = EL (eV)"""
    t = np.array(uc[key]['T'])[2:-2]; cb = np.array(uc[key]['Cb'])[2:-2]
    s = np.gradient(np.log(cb), np.log(t))
    return 1 + np.interp(np.log(KB*np.asarray(TK)/EL), np.log(t), s)
def n_window(key, EL, lo, hi):
    """exponent of a power-law fit rho - rho_0 = A T^n over [lo, hi] to the predicted T*C_b(T) (rho_0 = 0)"""
    T = np.linspace(lo, hi, 300); t = np.array(uc[key]['T'])[2:-2]; cb = np.array(uc[key]['Cb'])[2:-2]
    C = np.exp(np.interp(np.log(KB*T/EL), np.log(t), np.log(cb)))
    return np.polyfit(np.log(T), np.log(T*C), 1)[0]
say('(1) Predicted effective exponent n of rho - rho_0 (conformal bath: n = 1 exactly)')
WINDOWS = [(1.0, 26.7), (1.44, 30.0), (40.0, 85.0)]
pred = {}
for key, lab in CASES.items():
    for EL in (1.0, 10.0):
        ns = [n_window(key, EL, lo, hi) for lo, hi in WINDOWS]
        pred[f'{key}_EL{EL:g}'] = ns
        say(f'  {lab:32s} E_L/gamma = {EL:4.0f} eV: ' + '   '.join(f'{lo:g}-{hi:g} K: n = {n:.3f}' for (lo, hi), n in zip(WINDOWS, ns)))
allp = np.array(list(pred.values()))
say(f'  range over all cases: ' + '   '.join(f'{lo:g}-{hi:g} K: {allp[:, i].min():.3f}-{allp[:, i].max():.3f}' for i, (lo, hi) in enumerate(WINDOWS)))
say()

# ---------------------------------------------------------------- (2) data
plaw = lambda T, r0, A, n: r0 + A*T**n
def fit_n(T, y, lo, hi):
    m = (T >= lo) & (T <= hi)
    p, cv = curve_fit(plaw, T[m], y[m], p0=[y[m].min(), (y[m].max()-y[m].min())/(hi-lo), 1.0], maxfev=20000)
    return p[2], np.sqrt(cv[2, 2]), m.sum()
def local_slopes(T, y, edges):
    out = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        m = (T >= lo) & (T < hi)
        if m.sum() > 8:
            p, cv = np.polyfit(T[m], y[m], 1, cov=True); out.append((np.sqrt(lo*hi), p[0], np.sqrt(cv[0, 0])))
    return out

lg = json.load(open(D('legros2019_fig1_rho.json')))
mi = json.load(open(D('michon2023_fig2a_rho.json')))
summary = {'pred_windows': WINDOWS, 'pred': pred}

say('(2a) Nd-LSCO p = 0.24, H = 33 T (Daou 2009 via Legros 2019 Fig. 1c)')
T = np.array(lg['NdLSCO_p024_33T']['T']); y = np.array(lg['NdLSCO_p024_33T']['drho'])
ns = []
for lo in (1.0, 1.5, 2.0, 3.0):
    for hi in (20.0, 23.0, 26.7):
        n, dn, k = fit_n(T, y, lo, hi); ns.append(n)
        say(f'  window {lo:4.1f}-{hi:4.1f} K: n = {n:.4f} +- {dn:.4f} (fit error, {k} points)')
# block bootstrap over 1-K blocks for the main window
rng = np.random.default_rng(1); m = (T >= 1.0) & (T <= 26.7); Tm, ym = T[m], y[m]
blocks = [np.where((Tm >= a) & (Tm < a+1))[0] for a in np.arange(1.0, 26.7, 1.0)]
bs = []
for _ in range(300):
    idx = np.concatenate([blocks[i] for i in rng.integers(0, len(blocks), len(blocks))])
    try: bs.append(curve_fit(plaw, Tm[idx], ym[idx], p0=[0, 0.47, 1.0], maxfev=20000)[0][2])
    except RuntimeError: pass
nd_n = (float(min(ns)), float(max(ns))); nd_bs = float(np.std(bs))
say(f'  -> n = {nd_n[0]:.3f}-{nd_n[1]:.3f} across windows; block-bootstrap s.d. (1-K blocks, 1-26.7 K) = {nd_bs:.3f}')
for Tc, s, ds in local_slopes(T, y, [1, 3, 6, 10, 15, 20, 26.7]):
    say(f'     local slope at {Tc:5.1f} K: {s:.3f} +- {ds:.3f} uOhm cm/K')
summary['NdLSCO'] = dict(window=[1.0, 26.7], n=nd_n, n_bootstrap_sd=nd_bs)

say('(2b) LSCO p = 0.23, H = 48 T (Cooper 2009 via Legros 2019 Fig. 1b; 16 points)')
T = np.array(lg['LSCO_p023_48T']['T']); y = np.array(lg['LSCO_p023_48T']['drho'])
n, dn, k = fit_n(T, y, 1.0, 31.0)
say(f'  1.4-30 K: n = {n:.3f} +- {dn:.3f}')
summary['LSCO23'] = dict(window=[1.44, 30.05], n=float(n), dn=float(dn))
say(f'     vs prediction over 1.44-30 K (0.954-0.988; fixed amplitude 0.954-0.972): '
    f'{(n-0.988)/dn:.1f} s.d. above the top, {(n-0.972)/dn:.1f} above the fixed-amplitude top, {(n-0.954)/dn:.1f} above the bottom')

say('(2c) LSCO p = 0.24 (Michon 2023 Fig. 2a; the optics crystal)')
for key in ('rho_0T', 'rho_16T'):
    T = np.array(mi[key]['T']); y = np.array(mi[key]['rho']); u, ix = np.unique(T, return_index=True); T, y = T[ix], y[ix]
    edges = [20, 30, 40, 55, 70, 85, 100, 120] + ([150, 200, 263] if key == 'rho_0T' else [])
    say(f'  {key}: local slopes (uOhm cm/K): ' + ', '.join(f'{Tc:.0f} K {s:.4f}' for Tc, s, ds in local_slopes(T, y, edges)))
    for lo, hi in ((40, 85), (35, 75), (30, 80), (20, 80)):
        if key == 'rho_0T' and lo < 30: continue
        n, dn, k = fit_n(T, y, lo, hi); say(f'     {lo}-{hi} K: n = {n:.4f} +- {dn:.4f}')
        if (lo, hi) == (40, 85): summary[f'LSCO24_{key}'] = dict(window=[40, 85], n=float(n), dn=float(dn))
T = np.array(mi['rho_0T']['T']); y = np.array(mi['rho_0T']['rho']); u, ix = np.unique(T, return_index=True); T, y = T[ix], y[ix]
ns24 = [fit_n(T, y, lo, hi)[0] for lo, hi in ((40, 85), (35, 75), (30, 80), (45, 90), (40, 70))]
summary['LSCO24_0T_range'] = [float(min(ns24)), float(max(ns24))]
say(f'  0 T: n = {min(ns24):.3f}-{max(ns24):.3f} across windows 30-90 K')
# the predicted window itself, 40-85 K, and its sub-windows (Fig. 2b and the Letter)
ns4085 = [fit_n(T, y, lo, hi)[0] for lo, hi in ((40, 85), (40, 70), (40, 75), (40, 80), (45, 80), (45, 85), (50, 80), (50, 85))]
summary['LSCO24_0T_40_85_range'] = [float(min(ns4085)), float(max(ns4085))]
say(f'  0 T, predicted window 40-85 K: n = {ns4085[0]:.3f} (sub-windows {min(ns4085):.3f}-{max(ns4085):.3f}); prediction 0.934-0.986')
# model dependence of a fit that includes the superlinear term above 85 K
def c_fit(bg, lo, hi):
    m = (T >= lo) & (T <= hi)
    if bg == 'T2': f = lambda T, r0, A, c, A2: r0 + A*T*(1 + c*np.log10(55/T)) + A2*T**2
    else: f = lambda T, r0, A, c, A3: r0 + A*T*(1 + c*np.log10(55/T)) + A3*T**3
    p, cv = curve_fit(f, T[m], y[m], p0=[13, 0.63, 0, 1e-5], maxfev=20000); return p[2]
cs = [c_fit(bg, 40, hi) for bg in ('T2', 'T3') for hi in (150, 200, 263)]
say(f'  0 T, 40-150/200/263 K with a T^2 or T^3 background: fitted drift c = {min(cs):+.2f} to {max(cs):+.2f} per decade '
    f'(prediction +0.03 to +0.15): the answer is set by the background model')
summary['LSCO24_bg_c'] = [float(min(cs)), float(max(cs))]
say()

say('(3) Optics, LSCO p = 0.24 (Michon 2023 Fig. 2b): zero-frequency value per temperature')
oc = json.load(open(D('michon2023_fig2b_curves.json')))
for c in oc:
    a = np.unique(np.array(c['pts']), axis=0); x, yy = a[:, 0], a[:, 1]
    say(f"  {c['T']:3d} K: median hbar/(tau kT) over hbar w/kT < 1.5 = {np.median(yy[x < 1.5]):.2f}; over 3-6: {np.median(yy[(x >= 3) & (x < 6)]):.2f}")
say('  The 100-300 K values scatter by +-15% without a monotonic trend, and the 50-75 K curves fail to collapse below')
say('  hbar w/kT ~ 10; the predicted change of C_b over 100-300 K is 2-7%, so the optics cannot test the drift.')
json.dump(summary, open(D('drift_summary.json'), 'w'), indent=1)
open(os.path.join(H, '..', 'outputs', 'drift_test.txt'), 'w').write('\n'.join(lines) + '\n')
