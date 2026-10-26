"""Extract the curves of Fig. 2a (resistivity) and Fig. 2b (optical scattering rate) of
Michon et al., Nat. Commun. 14, 3033 (2023), arXiv:2205.04030, from the vector paths of the arXiv
source file fig02.pdf. Rebuilds data/michon2023_fig2b_curves.json (the original extraction code was lost, and the old file
carried a calibration offset: it set y = 0 at the centre of the '0' tick label, 2.2 pt below the axis,
so every hbar/(tau kT) value was 0.906 too high and every hbar w/kT value 0.018 too low; slopes are
unaffected) and adds data/michon2023_fig2a_rho.json (same crystal: dc resistivity at 0 T and 16 T).

Usage:  python3 extract_michon2023_fig2.py /path/to/fig02.pdf [outdir]
Get fig02.pdf from the arXiv source:  curl -L -o src.tgz https://arxiv.org/e-print/2205.04030 ; tar xzf src.tgz
Requires PyMuPDF (pip install pymupdf).

Axis calibration uses the frame lines and major tick marks of the PDF, not the tick-label text:
  panel a: x frame 37.21 (0 K) .. 207.28 (300 K); y frame 176.80 (0) .. 6.73 (0.2 mOhm cm)
  panel b: x frame 263.97 (0) .. 434.04 (100);    y frame 176.80 (0) .. 6.73 (70)
Check: the dashed guide line of panel a comes out as rho = 13.5 + 0.629 T (uOhm cm); the caption
quotes rho_0 = 12.2 uOhm cm and A = 0.63 uOhm cm/K, so the slope calibration is good to 0.2%
(the 1.3 uOhm cm offset is in the figure itself, about one point of the drawing).
"""
import sys, os, json
import numpy as np
import pymupdf

PANEL_A = dict(x0=37.21, x1=207.28, X0=0.0, X1=300.0, y0=176.80, y1=6.73, Y0=0.0, Y1=200.0)  # T in K, rho in uOhm cm
PANEL_B = dict(x0=263.97, x1=434.04, X0=0.0, X1=100.0, y0=176.80, y1=6.73, Y0=0.0, Y1=70.0)  # hbar w/kT, hbar/(tau kT)
# Fig. 2b legend: 40 K (green, leaves the frame), then 50 ... 300 K in this colour order
FIG2B_COLOURS = {50: (0.8438, 0.9219, 0.1562), 60: (1.0, 0.9258, 0.0), 75: (1.0, 0.7715, 0.0),
                 100: (1.0, 0.6172, 0.0), 150: (0.9648, 0.4648, 0.0), 200: (0.8105, 0.3105, 0.0),
                 250: (0.6562, 0.1562, 0.0), 300: (0.502, 0.002, 0.0)}


def calib(p, P):
    x = P['X0'] + (p[:, 0] - P['x0']) * (P['X1'] - P['X0']) / (P['x1'] - P['x0'])
    y = P['Y0'] + (p[:, 1] - P['y0']) * (P['Y1'] - P['Y0']) / (P['y1'] - P['y0'])
    return np.column_stack([x, y])


def path_points(d):
    pts = []
    for it in d['items']:
        if it[0] == 'l':
            pts += [(it[1].x, it[1].y), (it[2].x, it[2].y)]
    return np.array(pts)


def inside(p, P, pad=0.01):
    return ((p[:, 0] >= P['x0'] - pad) & (p[:, 0] <= P['x1'] + pad) &
            (p[:, 1] <= P['y0'] + pad) & (p[:, 1] >= P['y1'] - pad))


def same(c, ref, tol=2e-3):
    return c is not None and all(abs(a - b) < tol for a, b in zip(c, ref))


def main(pdf, outdir):
    drawings = pymupdf.open(pdf)[0].get_drawings()
    # ---- panel b: 8 curves, line width 0.75, one path per colour inside the frame
    curves = []
    for T, col in FIG2B_COLOURS.items():
        best = None
        for d in drawings:
            if same(d.get('color'), col) and abs(d['width'] - 0.75) < 0.01:
                p = path_points(d)
                m = inside(p, PANEL_B)
                if m.sum() > (0 if best is None else best[0]):
                    best = (m.sum(), p[m])
        a = calib(best[1], PANEL_B)
        curves.append(dict(T=T, color=list(col), n=int(len(a)), xmin=float(a[:, 0].min()),
                           xmax=float(a[:, 0].max()), pts=a.tolist()))
    # ---- panel a: zero field (black, width 2) and 16 T (red, width 1) inside the main frame
    rho = {}
    for key, col, w in (('rho_0T', (0, 0, 0), 2.0), ('rho_16T', (1, 0, 0), 1.0)):
        best = None
        for d in drawings:
            if same(d.get('color'), col) and abs(d['width'] - w) < 0.01 and len(d['items']) > 100:
                p = path_points(d)
                m = inside(p, PANEL_A)
                if m.sum() > 0 and p[:, 0].min() < PANEL_A['x0'] + 5 and (best is None or m.sum() > best[0]):
                    best = (m.sum(), p[m])
        a = calib(best[1], PANEL_A)
        a = a[np.argsort(a[:, 0])]
        rho[key] = dict(T=a[:, 0].tolist(), rho=a[:, 1].tolist())
    for d in drawings:   # dashed guide rho_0 + A T
        if same(d.get('color'), (0.5, 0.5, 0.5)) and len(d['items']) == 1 and d['rect'].y1 < 180:
            p = calib(path_points(d), PANEL_A)
            A = (p[1, 1] - p[0, 1]) / (p[1, 0] - p[0, 0])
            rho['guide_line'] = dict(rho0=float(p[0, 1] - A * p[0, 0]), A=float(A))
    rho['units'] = 'T in K, rho in micro-Ohm cm'
    os.makedirs(outdir, exist_ok=True)
    json.dump(curves, open(os.path.join(outdir, 'michon2023_fig2b_curves.json'), 'w'))
    json.dump(rho, open(os.path.join(outdir, 'michon2023_fig2a_rho.json'), 'w'))
    for c in curves:
        print(f"Fig 2b  T={c['T']:3d} K  n={c['n']}  x={c['xmin']:.3f}-{c['xmax']:.2f}")
    for k in ('rho_0T', 'rho_16T'):
        T = np.array(rho[k]['T']); print(f"Fig 2a  {k}: {len(T)} points, T = {T.min():.1f}-{T.max():.1f} K")
    print(f"guide line: rho = {rho['guide_line']['rho0']:.2f} + {rho['guide_line']['A']:.4f} T")


if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2] if len(sys.argv) > 2 else os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'data'))
