"""Every number quoted in the v6 Letter and Supplemental Material that is computed (not read from a paper),
with the section where it appears. Writes ../outputs/v6_numbers.txt. Run after drift_test_v6.py and band_dos.py.
Inputs read from published papers are listed in INPUTS below with their source; see ../README.md."""
import os, sys, json, math
import numpy as np
from scipy.special import lambertw

H = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, H)
from bath import w_conf, w_arc, S0, Cb, QS, CB_CONF, beta_of_Cb, r_of_Cb, Cb_of_beta
D = lambda f: os.path.join(H, '..', 'data', f)
O = lambda f: os.path.join(H, '..', 'outputs', f)
L = []
def say(s=''):
    print(s); L.append(s)

INPUTS = dict(
    gamma_b=3.89,         # mJ K^-2 mol^-1, Michon 2023 tight-binding band (t'/t = -0.17, m_b = 2.76 m_e)
    alpha_optics=5.0,     # hbar/(tau(0) kB T) quoted from the optics, Michon 2023 text
    alpha_dc=5.7,         # same from the measured A = 0.63 uOhm cm/K with Michon's conversion (5 <-> 0.55)
    lam_mass=0.46,        # twice Michon's g = 0.23 (optical mass logarithm)
    B=(1.9, 2.2, 2.5),    # Girod 2021, LSCO p = 0.24
    release_ratio=(1.47, 1.66),   # dc rate over the optical-mass-log scale, relative to conformal (read from Michon 2023)
    Q_Nd_B323={'in-plane': (7.41, 1.22), 'ADMR': (4.07, 1.45)},   # Nd-LSCO Q at B = 3.23 (v4 stress tests)
    B_Nd=2.87,            # Nd-LSCO logarithm above 2 K from the extracted curves (fit_michon2019_B.py)
    Q_Kondo={'Ce3Pd20Si6': (0.22, 0.32), 'CeCoIn5': (0.35, 0.35), 'YbRh2Si2': (0.43, 0.43), 'CeCu5.9Au0.1': (1.83, 1.83)},
)

say('== S2-S3. Bath functionals')
for name, w in (('tanh(x/2) conformal', w_conf), ('arctan, r=1.104', w_arc(1.104)), ('1-exp(-x)', lambda x: -np.expm1(-x)),
                ('tanh x', np.tanh), ('min(1,x)', lambda x: min(1.0, x)), ('arctan, r=0.45', w_arc(0.45)), ('arctan, r=0.30', w_arc(0.30))):
    b, c = S0(w), Cb(w); say(f'  {name:22s} beta = {b:.3f}  C_b = {c:.3f}  C_b/beta = {c/b:.3f}  R_opt = 2C_b/pi = {2*c/np.pi:.2f}')
fams = {'tanh(x/2s)': lambda s: (lambda x: np.tanh(x/(2*s))), '1-exp(-x/s)': lambda s: (lambda x: -np.expm1(-x/s)),
        'min(1,x/s)': lambda s: (lambda x: min(1.0, x/s)), 'x^2/(s^2+x^2)': lambda s: (lambda x: x*x/(s*s+x*x)),
        'tanh(x/2s)^2': lambda s: (lambda x: np.tanh(x/(2*s))**2), 'erf(x/s)': lambda s: (lambda x: math.erf(x/s))}
rat = []
for name, f in fams.items():
    for s in (0.5, 0.75, 1, 1.5, 2, 3):
        b = S0(f(s))
        if b >= 1.5: rat.append(Cb(f(s))/b)
ra = [Cb(w_arc(r))/S0(w_arc(r)) for r in np.linspace(0.3, 2, 12)]
say(f'  C_b/beta for beta >= 1.5: six families x six widths: {min(rat):.2f}-{max(rat):.2f} ({len(rat)} shapes); arctan r = 0.3-2: {min(ra):.2f}-{max(ra):.2f}')
say(f'  beta for arctan r = 0.30, 0.42, 0.45, 0.50, 0.53, 1, 2: ' + ', '.join(f'{S0(w_arc(r)):.2f}' for r in (0.30, 0.42, 0.45, 0.50, 0.53, 1, 2)))
say()
say('== S4. Conformal bath')
say(f'  C_b = 2 pi^3/7 zeta(3) = {CB_CONF:.4f};  Q_conf = 4 pi^3/21 zeta(3) = {QS:.4f};  I_w = pi^2/6 = {np.pi**2/6:.4f} -> b = 1 + 3 I_w/pi^2 = 1.5')
say()
say('== S5. Self-damped field: thermal-mass saddle (data/universal_curves.json)')
uc = json.load(open(D('universal_curves.json')))
allb, alld, allq, famb = [], [], [], []
for k, v in uc.items():
    T = np.array(v['T'])[2:-2]; m = (T >= 1e-5) & (T <= 1e-2)
    be = np.array(v['beta'])[2:-2][m]; de = np.array(v['delta'])[2:-2][m]; cb = np.array(v['Cb'])[2:-2][m]
    q = cb/(1+de); allb += list(be); alld += list(de); allq += list(q)
    if k.endswith('u1000'): famb += list(be)
    sl = np.gradient(np.log(np.array(v['Cb'])[2:-2]), np.log10(np.array(v['T'])[2:-2]))[m]
    db = np.gradient(np.array(v['beta'])[2:-2], np.log10(np.array(v['T'])[2:-2]))[m]
    say(f'  {k:12s}: beta {be.min():.2f}-{be.max():.2f}, delta {de.min():+.3f}..{de.max():+.3f}, eta = 1/(1+delta) {1/(1+de.max()):.3f}-{1/(1+de.min()):.3f}, '
        f'Q {q.min():.1f}-{q.max():.1f}, dbeta/dlog10T {db.min():+.2f}..{db.max():+.2f}, dlnC_b/dlog10T {sl.min():+.3f}..{sl.max():+.3f}')
say(f'  all cases, gamma T/E_L = 1e-5..1e-2: beta {min(allb):.2f}-{max(allb):.2f}, delta {min(alld):+.3f}..{max(alld):+.3f}, Q {min(allq):.1f}-{max(allq):.1f}')
say(f'  fixed amplitude (u -> inf): beta >= {min(famb):.2f}')
ell = 5.57; r = np.pi*lambertw(ell/np.pi).real/ell
say(f'  closed form r = pi W0(l/pi)/l: r({ell}) = {r:.3f}; a0/T = 2 pi exp(l - gamma_E) = {2*np.pi*np.exp(ell-np.euler_gamma):.0f}')
say(f'  identity: gamma_b lambda = (2 pi^2/3) g^2 N0^2 rho = pi rho gamma/3 with gamma = 2 pi g^2 N0^2 (checked symbolically in SM Sec. S5)')
say()
say('== S6. LSCO p = 0.24')
oc = json.load(open(D('michon2023_fig2b_curves.json')))
sl = []
for c in oc:
    p = np.array(c['pts']); m = (p[:, 0] > 20) & (p[:, 0] < min(80, p[:, 0].max()-1))
    if m.sum() > 50 and p[:, 0].max() > 45: sl.append(np.polyfit(p[m, 0], p[m, 1], 1)[0])
s = np.mean(sl); lam = 2*s/np.pi
say(f'  optical slope over hbar w/kT = 20-80 (50-100 K curves): {min(sl):.3f}-{max(sl):.3f}, mean {s:.3f} -> lambda = {lam:.3f}')
a0, a1 = INPUTS['alpha_optics'], INPUTS['alpha_dc']
Ropt = (a0/s, a1/s); cbs = (np.pi*Ropt[0]/2, np.pi*Ropt[1]/2)
say(f'  R_opt = {Ropt[0]:.2f}-{Ropt[1]:.2f} (conformal {2*CB_CONF/np.pi:.2f}; ratio {Ropt[0]/(2*CB_CONF/np.pi):.2f}-{Ropt[1]/(2*CB_CONF/np.pi):.2f}); '
    f'C_b = {cbs[0]:.1f}-{cbs[1]:.1f}; r = {r_of_Cb(cbs[1]):.2f}-{r_of_Cb(cbs[0]):.2f}; beta = {beta_of_Cb(cbs[0]):.2f}-{beta_of_Cb(cbs[1]):.2f}')
rel = [CB_CONF*x for x in INPUTS['release_ratio']]
say(f'  release route: C_b = {rel[0]:.1f}-{rel[1]:.1f}, r = {r_of_Cb(rel[1]):.2f}-{r_of_Cb(rel[0]):.2f}, beta = {beta_of_Cb(rel[0]):.2f}-{beta_of_Cb(rel[1]):.2f}')
say(f'  both routes: C_b = {rel[0]:.1f}-{cbs[1]:.1f}, beta = {beta_of_Cb(rel[0]):.2f}-{beta_of_Cb(cbs[1]):.2f}')
g = INPUTS['gamma_b']; Bs = INPUTS['B']
Q = [a*g/B for a in (a0, a1) for B in Bs]
say(f'  Q = alpha gamma_b/B: {a0*g/Bs[1]:.2f}-{a1*g/Bs[1]:.2f} at B = 2.2; {min(Q):.2f}-{max(Q):.2f} with B = 1.9-2.5')
for l in (lam, INPUTS['lam_mass']):
    say(f'  eta = gamma_b lambda/B with lambda = {l:.2f}: {g*l/Bs[1]:.3f} at B = 2.2; {g*l/Bs[2]:.3f}-{g*l/Bs[0]:.3f} with B = 1.9-2.5')
say(f'  check: eta = Q/C_b over Q 8.84-10.08 and C_b 10.8-14.3: {8.84/14.3:.2f}-{10.08/10.8:.2f}')
bd = open(O('band_dos.txt')).read().splitlines()
ref = None
for line in bd[1:]:
    toks = line.split()
    try: prodK = float(line[line.index(')')+1:].split()[2]); dos = float(line[line.index(')')+1:].split()[0])
    except (ValueError, IndexError): continue
    if ref is None: ref = (dos, prodK)
    say(f'  band {line[:38].strip():38s}: eta factor gamma_b-only {dos/ref[0]:.2f}, (m/m_e)K {prodK/ref[1]:.2f}')
say()
say('== S7. Nd-LSCO and other materials')
for key, (q, dq) in INPUTS['Q_Nd_B323'].items():
    qq, dd = q*3.23/INPUTS['B_Nd'], dq*3.23/INPUTS['B_Nd']
    say(f'  Nd-LSCO {key}: Q = {qq:.2f} +- {dd:.2f}; eta = Q/C_b for C_b 10.8-14.3: {qq/14.3:.2f}-{qq/10.8:.2f}; if conformal (C_b = 7.37): {qq/CB_CONF:.2f} +- {dd/CB_CONF:.2f}')
for k, (lo, hi) in INPUTS['Q_Kondo'].items():
    say(f'  {k}: eta = Q/C_b for C_b 7.37-12.5: {lo/12.5:.3f}-{hi/CB_CONF:.3f}')
say()
say('== S8. Drift test (drift_test_v6.py; data/drift_summary.json)')
sm = json.load(open(D('drift_summary.json')))
P = np.array(list(sm['pred'].values()))
say(f"  predicted n over 1-27 K: {P[:, 0].min():.3f}-{P[:, 0].max():.3f}; over 40-85 K: {P[:, 2].min():.3f}-{P[:, 2].max():.3f}; conformal: 1")
say(f"  Nd-LSCO p=0.24, 33 T: n = {sm['NdLSCO']['n'][0]:.3f}-{sm['NdLSCO']['n'][1]:.3f} (bootstrap s.d. {sm['NdLSCO']['n_bootstrap_sd']:.3f})")
say(f"  LSCO p=0.23, 48 T: n = {sm['LSCO23']['n']:.3f} +- {sm['LSCO23']['dn']:.3f}")
say(f"  LSCO p=0.24, 0 T, 30-90 K windows: n = {sm['LSCO24_0T_range'][0]:.3f}-{sm['LSCO24_0T_range'][1]:.3f}")
say(f"  LSCO p=0.24, 0 T, 40-263 K with a T^2/T^3 background: c = {sm['LSCO24_bg_c'][0]:+.2f}..{sm['LSCO24_bg_c'][1]:+.2f} per decade")
say()
say('== S7 and Fig. 3a. Eu-LSCO logarithm (data/michon2019_curves.json)')
c = json.load(open(D('michon2019_curves.json')))
T = np.array(c['Eu_p24']['T']); C = np.array(c['Eu_p24']['CT']); m = (T > 0.55) & (T < 10)
x = np.log(T[m]); A = np.vstack([-x, np.ones_like(x)]).T; (B, a0_), *_ = np.linalg.lstsq(A, C[m], rcond=None)
pp = np.polyfit(x, np.log(C[m]), 1)
say(f'  0.55-10 K: B = {B:.3f}, rms {np.std(C[m]-(a0_-B*x)):.3f}; power law exponent {pp[0]:.3f}, rms {np.std(C[m]-np.exp(np.polyval(pp, x))):.3f}')
open(O('v6_numbers.txt'), 'w').write('\n'.join(L) + '\n')
