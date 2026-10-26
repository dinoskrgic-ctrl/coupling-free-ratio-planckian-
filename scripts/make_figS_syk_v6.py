"""Supplemental Fig. S_syk: the two-band SYK check of Q = Q_conf (restyled repro/syk_two_band/final.py).
Reads ../syk_two_band/out/*.npy, written by ../syk_two_band/run_matsu*.py (Matsubara free energy) and
run_real.py / run_more.py (real-axis Kubo rate); see ../syk_two_band/run_syk.sh to recompute them.
Writes figS_syk.pdf/png next to the .tex files and ../outputs/figS_syk_table.txt."""
import numpy as np, os, sys
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
import scienceplots  # noqa: F401
plt.style.use(["science", "no-latex"]); plt.rcParams.update({"font.size": 8, "mathtext.fontset": "cm", "savefig.dpi": 300})
HERE = os.path.dirname(os.path.abspath(__file__))
SYK = os.path.join(HERE, '..', 'syk_two_band')
sys.path.insert(0, SYK)
from analyze import fitB, Cb, Qs
out = os.path.join(SYK, 'out') + '/'
D = 1.0; gb = 2 * np.pi / 3

def alpha0(T, a):
    m = T <= 1e-3
    vals = [np.polyfit(T[m], a[m], 1)[-1], np.polyfit(T[m], a[m], 2)[-1], np.polyfit(T, a, 3)[-1]]
    return np.mean(vals), np.ptp(vals)

def Bfit(F):
    Bs = [fitB(F[:, 0], F[:, 1], tm, True)[0] for tm in (1e-3, 2e-3, 4e-3)]
    return np.mean(Bs), np.ptp(Bs)

def local_lamC(F, Tgrid):
    """-(2/3gb) d(C/T)/dlnT from the fitted F = a + bT^2 + cT^2 lnT + dT^3 + eT^3 lnT."""
    _, _, c = fitB(F[:, 0], F[:, 1], 4e-3, True)
    a, b, cc, d, e = c
    # C/T = -F''(T):  F'' = 2b + cc(2lnT+3) + 6dT + e(6T lnT + 5T)
    CT = -(2 * b + cc * (2 * np.log(Tgrid) + 3) + 6 * d * Tgrid + e * (6 * Tgrid * np.log(Tgrid) + 5 * Tgrid))
    dCT = -(2 * cc + 6 * d * Tgrid + e * (6 * Tgrid * np.log(Tgrid) + 11 * Tgrid))   # d(C/T)/dlnT
    return CT, -dCT * 2 / (3 * gb)

rows = []
d = np.load(out + 'real_alpha.npy'); T = d[:, 0]
cases = [(0.17, 1, 'F_U0.17_dt0.02.npy'), (0.3, 4, 'F_U0.3_dt0.02.npy'), (0.5, 7, 'F_U0.5_dt0.02.npy')]
fig, ax = plt.subplots(1, 2, figsize=(6.8, 2.5))
cols = ['#1f4e79', '#2b8cbe', '#a6bddb']
for (U, j, fF), col in zip(cases, cols):
    F = np.load(out + fF); a0, da = alpha0(T, d[:, j]); B, dB = Bfit(F)
    Q = a0 * gb / B
    rows.append(('frozen bath', 0.05, U, 0, 0.02, a0 / Cb, 2 * B / (3 * gb), Q, Q / Qs))
    Tg = np.exp(np.linspace(np.log(1e-4), np.log(4e-3), 100))
    CT, lamC = local_lamC(F, Tg)
    ax[0].plot(Tg, lamC, color=col, lw=1.2, label=f'$U={U}$')
    ax[0].plot(T, d[:, j] / Cb, 'o', color=col, ms=3)
    ax[0].plot([7e-5], [a0 / Cb], '*', color=col, ms=7, mec='k', mew=0.4)
    ax[0].plot([7e-5], [2 * B / (3 * gb)], '_', color='k', ms=9, mew=1.0)
    ax[1].bar(len(rows) - 1, Q / Qs, color=col, width=0.6)
extra = [('frozen bath, $U_f=0.1$', 'F_U0.3_dt0.02_r0.0_Uf0.1.npy', 'real_alpha_uf0.1.npy', 0.1, 0.0, 0.02),
         ('frozen bath, $d\\tau=0.01$', 'F_U0.3_dt0.01_r0.0_Uf0.05.npy', None, 0.05, 0.0, 0.01),
         ('feedback, $N_c=N_f$', 'F_U0.3_dt0.02_r1.0_Uf0.05.npy', 'real_alpha_fc.npy', 0.05, 1.0, 0.02)]
for lab, fF, fa, Uf, r, dt in extra:
    if not os.path.exists(out + fF): continue
    F = np.load(out + fF); B, dB = Bfit(F)
    if fa is None: a = d[:, 4]; Ta = T
    else:
        if not os.path.exists(out + fa): continue
        dd = np.load(out + fa); Ta, a = dd[:, 0], dd[:, 1]
    if len(Ta) < 8: continue
    a0, da = alpha0(Ta, a); Q = a0 * gb / B
    rows.append((lab, Uf, 0.3, r, dt, a0 / Cb, 2 * B / (3 * gb), Q, Q / Qs))
    ax[1].bar(len(rows) - 1, Q / Qs, color='0.6', width=0.6)
ax[0].set_xscale('log'); ax[0].set_xlabel('$T/D$  (stars, bars: $T\\to0$)'); ax[0].set_ylabel('$\\lambda$')
ax[0].text(1.1e-4, 1.45, 'lines: $\\frac{2}{3\\gamma_b}\\,|d(C/T)/d\\ln T|$\ndots: $\\hbar/(\\tau_{\\rm tr}k_BT\\,C_b)$', fontsize=6.5)
ax[0].legend(fontsize=6.5, frameon=False, loc='center right')
ax[1].axhline(1, color='k', lw=0.8); ax[1].set_ylim(0.99, 1.01); ax[1].set_yticks([0.99, 0.995, 1.0, 1.005, 1.01])
ax[1].set_xticks(range(len(rows))); labs = []
for i, r in enumerate(rows):
    tag = '' if i < 3 else {0.1: '$U_f=0.1$', }.get(r[1], '') or ('$d\\tau/2$' if r[4] == 0.01 else '$N_c=N_f$')
    labs.append((tag + '\n' if tag else '') + f'$\\lambda={r[5]:.2f}$')
ax[1].set_xticklabels(labs, fontsize=6)
ax[1].set_ylabel(r'$Q/Q_{\rm conf}$')
fig.tight_layout(pad=0.3)
fig.savefig(os.path.join(HERE, '..', '..', 'figS_syk.pdf')); fig.savefig(os.path.join(HERE, '..', '..', 'figS_syk.png'))
with open(os.path.join(HERE, '..', 'outputs', 'figS_syk_table.txt'), 'w') as fh:
    fh.write('case | Uf/D | U | r=Nc/Nf | dtau | lam_tr=alpha/C_b | lam_C=2B/(3 gamma_b) | Q | Q/Q*\n')
    for r in rows:
        fh.write(f"{r[0]} | {r[1]} | {r[2]} | {r[3]} | {r[4]} | {r[5]:.4f} | {r[6]:.4f} | {r[7]:.4f} | {r[8]:.4f}\n")
print(open(os.path.join(HERE, '..', 'outputs', 'figS_syk_table.txt')).read())
