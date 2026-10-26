"""Figures for the v6 Letter (entropy_prl_v6.tex) and its Supplemental Material.
Single-column (3.375 in) figures in the SciencePlots 'science' style, with no titles inside the figures.

  fig1_map       the (beta, eta) map with contours of Q = C_b(beta) eta
  fig2_field     (a) beta of the self-damped field from the thermal-mass saddle; (b) the drift test: effective resistivity
                 exponent predicted by the field against three published curves (conformal bath: exactly 1)
  fig3_data      (a) Eu-/Nd-LSCO C_el/T at p = 0.24 with log and power-law fits; (b) LSCO p = 0.24 optical rate against
                 the local-bath model
  figS_spectra, figS_exponents, figS_materials, figS_drift   (Supplemental Material)
  figS_syk comes from make_figS_syk_v6.py.

Inputs: ../data/*.json (see ../README.md for provenance). Outputs: PDF and PNG in ../../ (next to the .tex files).
Needs numpy, scipy, matplotlib and SciencePlots (pip install SciencePlots).
STYLE uses 'no-latex' so that it runs without a TeX installation; with LaTeX installed, set STYLE = ["science"]."""
import json, os, sys
import numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import scienceplots  # noqa: F401  (registers the "science" styles)
from scipy.integrate import quad

H = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, H)
from bath import w_conf, w_arc, S0, Cb, optical_model, QS, CB_CONF, arctan_table, beta_of_Cb

D = lambda f: os.path.join(H, "..", "data", f)
OUT = lambda f: os.path.join(H, "..", "..", f)
STYLE = ["science", "no-latex"]
plt.style.use(STYLE)
plt.rcParams.update({"font.size": 8, "axes.labelsize": 8, "xtick.labelsize": 7, "ytick.labelsize": 7,
                     "legend.fontsize": 6.5, "mathtext.fontset": "cm", "savefig.dpi": 300})
CONF, BOS, KON, DATA, ND = "#009E73", "#0072B2", "#D55E00", "0.2", "#CC79A7"
W1 = 3.375
KB = 8.617e-5

# ---- inputs shared with v6_numbers.py
GAMMA_B = 3.89                 # LSCO p = 0.24 band Sommerfeld coefficient, mJ K^-2 mol^-1 (Michon 2023 band, m_b = 2.76 m_e)
ALPHA = (5.0, 5.7)             # hbar/(tau kB T) at omega -> 0: optics (5) and dc resistivity (5.7)
LAM = (0.40, 0.46)             # optical slope (2 x 0.625/pi) and twice Michon's g = 0.23
B_LSCO = (1.9, 2.2, 2.5)       # Girod 2021, B = 2.2 +- 0.3
Q_ND = {"in-plane": 8.34, "ADMR": 4.58}   # Nd-LSCO p = 0.24 with B = 2.87 (SM Sec. S6)
CB_LSCO = (10.8, 14.3)         # release route (low) to slope route with alpha = 5.7 (high)


def save(fig, name):
    fig.savefig(OUT(name + ".pdf")); fig.savefig(OUT(name + ".png")); plt.close(fig)


def label(ax, s, x=0.03, y=0.95, **kw):
    ax.text(x, y, s, transform=ax.transAxes, fontsize=8, va="top", ha="left", **kw)


def saddle():
    return json.load(open(D("universal_curves.json")))


def lsco_samples(n=25):
    """(beta, eta) for alpha, lambda, B spanning their ranges: beta from C_b = alpha/lambda, eta = gamma_b lambda / B"""
    out = []
    for a in np.linspace(*ALPHA, n):
        for l in np.linspace(*LAM, n):
            for B in np.linspace(B_LSCO[0], B_LSCO[2], 7):
                out.append((float(beta_of_Cb(a/l)), GAMMA_B*l/B, B))
    return np.array(out)


# ------------------------------------------------------------------ Fig. 1: the (beta, eta) map
def fig1():
    rs, be, cb = arctan_table()
    o = np.argsort(be); be, cb = be[o], cb[o]
    fig, a = plt.subplots(figsize=(W1, 2.9))
    xb = np.linspace(1.4, 5.25, 220); cbx = np.interp(xb, be, cb)
    a.axhspan(1.0, 1.3, color="0.88", lw=0, zorder=0)
    a.text(5.2, 1.2, r"$\eta>1$: no local bath", fontsize=6.3, color="0.35", va="center", ha="right")
    for Q, lab in ((2, "2"), (QS, "4.9"), (7, "7"), (10, "10"), (13, "13")):
        y = Q/cbx; m = y < 1.27
        a.plot(xb[m], y[m], color="0.6", lw=0.55, ls="--", zorder=1)
        if Q < 5:   # label at the left edge
            xl = 1.47; a.text(xl, Q/np.interp(xl, xb, cbx)+0.012, rf"$Q={lab}$", fontsize=5.6, color="0.4", ha="left", va="bottom")
        else:       # label where the contour leaves the top of the frame
            xl = float(np.interp(1.22, (Q/cbx)[::-1], xb[::-1]))
            a.text(xl+0.04, 1.22, rf"$Q={lab}$", fontsize=5.6, color="0.4", ha="left", va="center")
    # Kondo lattices: eta only (beta unmeasured)
    a.axhspan(0.02, 0.25, color=KON, alpha=0.10, lw=0, zorder=0)
    a.text(1.47, 0.135, r"Kondo lattices ($\beta$ unmeasured)", fontsize=6.2, color=KON, va="center")
    # photoemission range of beta in Bi-based cuprates
    a.plot([3.0, 4.0], [0.282, 0.282], color="0.45", lw=2.2, solid_capstyle="butt", zorder=2)
    a.text(4.08, 0.282, r"photoemission $\beta$ (Bi-based)", fontsize=5.6, color="0.35", ha="left", va="center",
           bbox=dict(fc="white", ec="none", pad=0.4), zorder=6)
    # self-damped field: (beta, 1/(1+delta)) along the saddle
    d = saddle()
    for k in d:
        T = np.array(d[k]["T"])[2:-2]; m = (T >= 1e-5) & (T <= 1e-2)
        b_ = np.array(d[k]["beta"])[2:-2][m]; e_ = 1/(1+np.array(d[k]["delta"])[2:-2][m])
        a.plot(b_, e_, color=BOS, lw=2.0 if k == "Ec1.0_u1000" else 0.9, alpha=0.9, solid_capstyle="round", zorder=4)
    a.text(2.38, 1.045, "self-damped 2D field", color=BOS, fontsize=6.3, ha="left", va="bottom")
    # LSCO p = 0.24
    s = lsco_samples()
    from scipy.spatial import ConvexHull
    hull = ConvexHull(s[:, :2]); a.fill(*s[hull.vertices, :2].T, fc=DATA, ec="none", alpha=0.13, zorder=2)
    c = s[np.isclose(s[:, 2], 2.2)]; hc = ConvexHull(c[:, :2])
    a.fill(*c[hc.vertices, :2].T, fc="none", ec=DATA, lw=1.0, zorder=5)
    a.text(4.17, 0.84, "LSCO\n$p=0.24$", color=DATA, fontsize=6.2, ha="left", va="center",
           bbox=dict(fc="white", ec="none", pad=0.4), zorder=6)
    # Nd-LSCO p = 0.24, two transport routes, beta assumed equal to LSCO's
    for key, q in Q_ND.items():
        lo, hi = q/CB_LSCO[1], q/CB_LSCO[0]
        a.add_patch(plt.Rectangle((3.0, lo), 1.1, hi-lo, fc="none", ec=ND, lw=0.8, ls=(0, (2.5, 1.5)), zorder=5))
        a.text(4.17, (lo+hi)/2, f"Nd-LSCO, {key}", color=ND, fontsize=5.8, ha="left", va="center",
               bbox=dict(fc="white", ec="none", pad=0.4), zorder=6)
    a.plot(2.0, 2/3, "*", color=CONF, ms=10, mec="k", mew=0.4, zorder=7)
    a.text(1.47, 0.70, "conformal bath\n" r"$\beta=2$, $\eta=2/3$, $Q=4.91$", color=CONF, fontsize=6.2, ha="left", va="bottom")
    a.set_xlim(1.4, 5.25); a.set_ylim(0, 1.27)
    a.set_xlabel(r"exponent $\beta$ (photoemission, or optics via $R_{\rm opt}=2C_b/\pi$)")
    a.set_ylabel(r"entropy share $\eta=\gamma_b\lambda/B$")
    save(fig, "fig1_map")
    return s


# ------------------------------------------------------------------ Fig. 2: the self-damped field and the drift test
def fig2():
    d = saddle(); sm = json.load(open(D("drift_summary.json")))
    fig, ax = plt.subplots(2, 1, figsize=(W1, 4.1))
    b = ax[0]
    for key, lab, ls, lw in (("Ec1.0_u0.3", r"$u\rho=0.3$", ":", 1.1), ("Ec1.0_u1", r"$u\rho=1$", "--", 1.1),
                             ("Ec1.0_u1000", r"$u\rho\to\infty$ (fixed amplitude)", "-", 1.7)):
        T = np.array(d[key]["T"])[2:-2]; b.plot(T, np.array(d[key]["beta"])[2:-2], color=BOS, ls=ls, lw=lw, label=lab)
    for key in ("Ec0.3_u1000", "Ec3.0_u1000"):
        T = np.array(d[key]["T"])[2:-2]; b.plot(T, np.array(d[key]["beta"])[2:-2], color=BOS, lw=0.5, alpha=0.6)
    b.axhspan(3.0, 4.0, color="0.5", alpha=0.15, lw=0)
    b.axhline(2.0, color=CONF, lw=1.4)
    b.text(1.5e-2, 2.07, "conformal", color=CONF, fontsize=6.3, va="bottom", ha="right")
    b.set_xscale("log"); b.set_ylim(1.6, 5); b.set_xlabel(r"$\gamma k_BT/E_L$"); b.set_ylabel(r"$\beta$")
    b.legend(frameon=False, loc="lower left", bbox_to_anchor=(0.0, 0.1)); label(b, "(a)", x=0.03, y=0.95)
    e = ax[1]
    TK = np.logspace(0, np.log10(300), 200)
    t = None
    for key in ("Ec1.0_u0.3", "Ec1.0_u1", "Ec1.0_u1000", "Ec0.3_u1000", "Ec3.0_u1000"):
        tt = np.array(d[key]["T"])[2:-2]; cb = np.array(d[key]["Cb"])[2:-2]; sl = np.gradient(np.log(cb), np.log(tt))
        curves = [1 + np.interp(np.log(KB*TK/EL), np.log(tt), sl) for EL in (1.0, 10.0)]
        lo, hi = np.minimum(*curves), np.maximum(*curves)
        t = (lo, hi) if t is None else (np.minimum(t[0], lo), np.maximum(t[1], hi))
    e.fill_between(TK, t[0], t[1], color=BOS, alpha=0.25, lw=0, label="self-damped field (prediction)")
    e.axhline(1.0, color=CONF, lw=1.4, label="conformal bath")
    rows = [("Nd-LSCO $p=0.24$, 33 T", sm["NdLSCO"]["window"], sm["NdLSCO"]["n"], DATA),
            ("LSCO $p=0.23$, 48 T", sm["LSCO23"]["window"], (sm["LSCO23"]["n"]-sm["LSCO23"]["dn"], sm["LSCO23"]["n"]+sm["LSCO23"]["dn"]), "0.45"),
            ("LSCO $p=0.24$, 0 T", [40, 85], sm["LSCO24_0T_40_85_range"], KON)]
    for lab, (lo, hi), (n0, n1), col in rows:
        e.add_patch(plt.Rectangle((lo, n0), hi-lo, n1-n0, fc=col, ec=col, alpha=0.35, lw=0.6, zorder=3))
        e.plot([], [], "s", color=col, alpha=0.5, ms=5, label=lab)
    e.set_xscale("log"); e.set_xlim(1, 300); e.set_ylim(0.90, 1.12)
    e.set_xlabel(r"$T$ (K)"); e.set_ylabel(r"$n=d\ln(\rho-\rho_0)/d\ln T$")
    e.legend(frameon=False, loc="upper right", ncol=2, columnspacing=0.8, handlelength=1.4); label(e, "(b)", x=0.03, y=0.12)
    fig.tight_layout(pad=0.3, h_pad=0.6); save(fig, "fig2_field")


# ------------------------------------------------------------------ Fig. 3: published data at p = 0.24
def fig3():
    fig, ax = plt.subplots(2, 1, figsize=(W1, 4.3))
    c = json.load(open(D("michon2019_curves.json")))
    a = ax[0]
    for key, col, lab in (("Nd_p24", "0.6", "Nd-LSCO"), ("Eu_p24", DATA, "Eu-LSCO")):
        T = np.array(c[key]["T"]); C = np.array(c[key]["CT"]); o = np.argsort(T); T, C = T[o], C[o]
        keep = (T > 0.5) & (T < 10.2); a.plot(T[keep][::6], C[keep][::6], ".", ms=1.5, color=col, label=f"{lab}, $p=0.24$")
    T = np.array(c["Eu_p24"]["T"]); C = np.array(c["Eu_p24"]["CT"]); m = (T > 0.55) & (T < 10)
    x = np.log(T[m]); A = np.vstack([-x, np.ones_like(x)]).T; (B, a0), *_ = np.linalg.lstsq(A, C[m], rcond=None)
    pp = np.polyfit(x, np.log(C[m]), 1); tt = np.logspace(np.log10(0.5), 1, 100)
    a.plot(tt, a0 - B*np.log(tt), color=BOS, lw=1.2, label=rf"$B\ln(T_0/T)$, $B={B:.2f}$")
    a.plot(tt, np.exp(np.polyval(pp, np.log(tt))), color=KON, lw=1.0, ls="--", label=rf"$\propto T^{{{pp[0]:.2f}}}$")
    a.set_xscale("log"); a.set_xlabel(r"$T$ (K)"); a.set_ylabel(r"$C_{\rm el}/T$ (mJ K$^{-2}$ mol$^{-1}$)")
    a.set_xticks([0.5, 1, 2, 5, 10]); a.set_xticklabels(["0.5", "1", "2", "5", "10"])
    a.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter()); a.set_ylim(10, 23.5)
    a.legend(frameon=False, loc="upper right"); label(a, "(a)", x=0.03, y=0.12)
    b = ax[1]; dd = json.load(open(D("michon2023_fig2b_curves.json")))
    slopes = []
    for k, cc in enumerate(dd):
        p = np.array(cc["pts"]); p = p[np.argsort(p[:, 0], kind="stable")]
        b.plot(p[::5, 0], p[::5, 1], color=plt.cm.Greys(0.35+0.08*k), lw=0.7, label="LSCO $p=0.24$, 50–300 K" if k == 0 else None)
        mm = (p[:, 0] > 20) & (p[:, 0] < min(80, p[:, 0].max()-1))
        if mm.sum() > 50 and p[:, 0].max() > 45: slopes.append(np.polyfit(p[mm, 0], p[mm, 1], 1)[0])
    s = float(np.mean(slopes)); lam = 2*s/np.pi
    xs = np.array([0.3, 1, 2, 4, 7, 10, 15, 20, 30, 40, 60, 80])
    for w, col, lab, ls in ((w_conf, CONF, "conformal", "-"), (w_arc(0.45), BOS, r"self-damped, $r=0.45$", "-"),
                            (w_arc(0.30), BOS, r"self-damped, $r=0.30$", "--")):
        b.plot(xs, optical_model(w, lam, xs), color=col, lw=1.2, ls=ls, label=lab)
    b.plot([0, 0], [5.0, 5.7], color=DATA, lw=4, solid_capstyle="butt", clip_on=False, zorder=5)
    b.set_xlim(0, 90); b.set_ylim(0, 62); b.set_xlabel(r"$\hbar\omega/k_BT$"); b.set_ylabel(r"$\hbar/\tau(\omega)k_BT$")
    b.legend(frameon=False, loc="upper left", bbox_to_anchor=(0.06, 1.0)); label(b, "(b)", x=0.88, y=0.12)
    fig.tight_layout(pad=0.3, h_pad=0.6); save(fig, "fig3_data")
    return B, pp[0], s, lam


# ------------------------------------------------------------------ Supplemental figures
def figS_spectra():
    fig, ax = plt.subplots(1, 2, figsize=(6.8, 2.3))
    a = ax[0]; x = np.linspace(0, 8, 400)
    a.plot(x, w_conf(x), color=CONF, lw=1.4, label=r"conformal, $\tanh(\Omega/2T)$")
    a.fill_between(x, w_arc(0.35)(x), w_arc(0.6)(x), color=BOS, alpha=0.2, lw=0)
    a.plot(x, w_arc(0.45)(x), color=BOS, lw=1.4, label=r"self-damped, $\frac{2}{\pi}\arctan(\Omega/rT)$, $r=0.45$")
    for xh, cc in ((2*np.arctanh(0.5), CONF), (0.45, BOS)):
        a.plot([xh, xh], [0, 0.5], ls=":", color=cc, lw=1); a.plot(xh, 0.5, "o", color=cc, ms=3)
    a.set_xlabel(r"$\hbar\Omega/k_BT$"); a.set_ylabel(r"$w=2\alpha^2F/\lambda$"); a.set_xlim(0, 8); a.set_ylim(0, 1.08)
    a.legend(loc="lower right", frameon=False); label(a, "(a)")
    c = ax[1]; x = np.linspace(0.01, 6, 400)
    c.plot(x, w_conf(x)/x, color=CONF, lw=1.4, label="conformal")
    c.plot(x, w_arc(0.45)(x)/x, color=BOS, lw=1.4, label=r"self-damped, $r=0.45$")
    c.fill_between(x, w_arc(0.35)(x)/x, w_arc(0.6)(x)/x, color=BOS, alpha=0.2, lw=0)
    c.set_ylim(0, 1.8); c.set_xlim(0, 6); c.set_xlabel(r"$\hbar\Omega/k_BT$"); c.set_ylabel(r"$T\chi''_{\rm loc}(\Omega)/\Omega$ (normalized)")
    c.legend(frameon=False, loc="upper right"); label(c, "(b)", x=0.85, y=0.12)
    fig.tight_layout(pad=0.3); save(fig, "figS_spectra")


def figS_materials():
    rows = [("two-band SYK (numerics)", QS*1.002, QS*1.000, QS*1.004, CONF),
            ("LSCO $p=0.24$", 9.4, 7.8, 11.7, DATA),
            ("Nd-LSCO $p=0.24$, in-plane", 8.34, 6.97, 9.71, ND),
            ("Nd-LSCO $p=0.24$, ADMR", 4.58, 2.95, 6.21, ND),
            (r"Sr$_3$Ru$_2$O$_7$ (7.9 T)", 4.7, 3.7, 6.1, "0.45"),
            (r"CeCu$_{5.9}$Au$_{0.1}$", 1.83, 0.39, 3.6, KON),
            (r"YbRh$_2$Si$_2$", 0.43, 0.17, 1.09, KON),
            (r"CeCoIn$_5$", 0.35, 0.18, 1.17, KON),
            (r"Ce$_3$Pd$_{20}$Si$_6$", 0.27, 0.22, 0.32, KON)]
    fig, a = plt.subplots(figsize=(W1, 2.4))
    a.axvspan(9, 13, color=BOS, alpha=0.12, lw=0); a.axvline(QS, color=CONF, lw=1.2)
    ys = []
    for i, (lab, q, lo, hi, cc) in enumerate(rows):
        y = len(rows)-1-i; ys.append(y); a.errorbar(q, y, xerr=[[q-lo], [hi-q]], fmt="o", color=cc, ms=3, capsize=1.5, lw=0.9)
    a.set_yticks(ys); a.set_yticklabels([r[0] for r in rows], fontsize=6); a.tick_params(axis="y", which="both", length=0)
    a.set_xscale("log"); a.set_xlim(0.1, 33); a.set_ylim(-0.7, len(rows)-0.3); a.set_xlabel(r"$Q=\alpha\gamma_b/B$")
    a.xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: f"{v:g}"))
    fig.tight_layout(pad=0.3); save(fig, "figS_materials")


def figS_exponents():
    fig, ax = plt.subplots(1, 2, figsize=(6.8, 2.4))
    a = ax[0]; s = np.linspace(-0.6, 0.6, 3)
    a.plot(1+s, s, color="k", lw=0.8); a.axvline(1.0, color=BOS, lw=5, alpha=0.25)
    pts = [("fractal, $d_s=4/3$", -1/3), ("$d=1.8$", -0.1), ("$d=2.2$", 0.1), ("membrane", 0.26),
           (r"long-range, $\sigma=1.5$", 0.33), ("3D, $z=2$", 0.5), ("2D, $z=2$", 0.0)]
    for lab, v in pts:
        a.plot(1+v, v, "o", color=KON if v != 0 else BOS, ms=3); a.text(1+v+0.025, v-0.05, lab, fontsize=5.8)
    a.set_xlabel(r"resistivity exponent $1+s$"); a.set_ylabel(r"$s$ in $C/T\propto T^{s}$"); a.set_xlim(0.55, 1.8)
    label(a, "(a)")
    b = ax[1]
    def rate(T, wmin, wmax=1.0):
        g = lambda W: 0.5*(2/np.expm1(W/T)+2/(np.exp(W/T)+1))
        return np.pi*quad(g, wmin, wmax, limit=400, points=[T])[0]
    Ts = np.logspace(-3, -1.3, 30)
    for wmin, ls in ((1e-4, "-"), (1e-3, "--")):
        r = np.array([rate(tt, wmin) for tt in Ts])
        b.plot(Ts, np.gradient(np.log(r), np.log(Ts)), color=KON, ls=ls, lw=1.1, label=rf"harmonic flat bath, $\omega_{{\min}}={wmin:g}$")
    b.axhline(1.0, color=BOS, lw=1.2, label=r"self-scaling bath ($\Gamma\propto T$)")
    b.set_xscale("log"); b.set_xlabel(r"$T/\omega_{\max}$"); b.set_ylabel(r"$d\ln\rho/d\ln T$"); b.set_ylim(0.9, 2.0)
    b.legend(frameon=False, loc="upper right"); label(b, "(b)")
    fig.tight_layout(pad=0.3); save(fig, "figS_exponents")


def figS_drift():
    lg = json.load(open(D("legros2019_fig1_rho.json"))); mi = json.load(open(D("michon2023_fig2a_rho.json")))
    def slopes(T, y, edges):
        out = []
        for lo, hi in zip(edges[:-1], edges[1:]):
            m = (T >= lo) & (T < hi)
            if m.sum() > 8:
                p, cv = np.polyfit(T[m], y[m], 1, cov=True); out.append((np.sqrt(lo*hi), p[0], np.sqrt(cv[0, 0])))
        return np.array(out)
    fig, ax = plt.subplots(1, 2, figsize=(6.8, 2.4))
    a = ax[0]
    T = np.array(lg["NdLSCO_p024_33T"]["T"]); y = np.array(lg["NdLSCO_p024_33T"]["drho"])
    s = slopes(T, y, [1, 2, 3, 4.5, 6, 8, 10, 13, 16, 20, 23, 26.7])
    a.errorbar(s[:, 0], s[:, 1], yerr=s[:, 2], fmt="o", ms=3, color=DATA, capsize=1.5, lw=0.8, label="Nd-LSCO $p=0.24$, 33 T")
    T2 = np.array(lg["LSCO_p023_48T"]["T"]); y2 = np.array(lg["LSCO_p023_48T"]["drho"])
    s2 = np.array([(np.sqrt(T2[i]*T2[i+3]), np.polyfit(T2[i:i+4], y2[i:i+4], 1)[0]) for i in range(0, len(T2)-3, 3)])
    a.plot(s2[:, 0], s2[:, 1]*0.465/np.mean(s2[:, 1]), "s", ms=3, color="0.55", label=r"LSCO $p=0.23$, 48 T (scaled)")
    a.set_xscale("log"); a.set_xlabel(r"$T$ (K)"); a.set_ylabel(r"$d\rho/dT$ ($\mu\Omega$ cm/K)"); a.set_ylim(0.33, 0.6)
    a.set_xticks([1, 2, 5, 10, 20]); a.set_xticklabels(["1", "2", "5", "10", "20"]); a.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    a.legend(frameon=False, loc="lower right"); label(a, "(a)")
    b = ax[1]
    for key, col, lab in (("rho_0T", KON, "0 T"), ("rho_16T", BOS, "16 T")):
        T = np.array(mi[key]["T"]); y = np.array(mi[key]["rho"]); u, ix = np.unique(T, return_index=True); T, y = T[ix], y[ix]
        edges = [25, 32, 40, 50, 60, 70, 85, 100, 120] + ([150, 200, 263] if key == "rho_0T" else [])
        s = slopes(T, y, edges); k = s[:, 0] > 30
        b.errorbar(s[k, 0], s[k, 1], yerr=s[k, 2], fmt="o-", ms=3, lw=0.8, color=col, capsize=1.5, label=f"LSCO $p=0.24$, {lab}")
    b.set_xscale("log"); b.set_xlabel(r"$T$ (K)"); b.set_ylabel(r"$d\rho/dT$ ($\mu\Omega$ cm/K)")
    b.set_xticks([30, 50, 100, 200]); b.set_xticklabels(["30", "50", "100", "200"]); b.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    b.legend(frameon=False, loc="upper left", bbox_to_anchor=(0.08, 1.0)); label(b, "(b)")
    fig.tight_layout(pad=0.3); save(fig, "figS_drift")


if __name__ == "__main__":
    s = fig1(); print(f"fig1: LSCO beta {s[:, 0].min():.2f}-{s[:, 0].max():.2f}, eta {s[:, 1].min():.2f}-{s[:, 1].max():.2f}")
    fig2(); print("fig2 written")
    B, sexp, slope, lam = fig3(); print(f"fig3: Eu-LSCO B={B:.3f}, power-law exponent {sexp:.3f}; optics slope {slope:.3f}, lambda {lam:.3f}")
    figS_spectra(); figS_materials(); figS_exponents(); figS_drift(); print("all v6 figures written")
