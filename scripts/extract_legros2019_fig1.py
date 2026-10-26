"""Extract two low-temperature resistivity curves from Fig. 1 of Legros et al., Nat. Phys. 15, 142 (2019),
arXiv:1805.02512 (the arXiv PDF has vector figures; Fig. 1 is on page 11, page index 10):
  panel c: Nd-LSCO, p = 0.24, H = 33 T, rho(T) - rho_0 (blue line; the data of Daou et al., Nat. Phys. 5, 31 (2009));
  panel b: LSCO, p = 0.23, H = 48 T, rho(T) - rho_0 (orange dots; the data of Cooper et al., Science 323, 603 (2009)).
Writes data/legros2019_fig1_rho.json.

Usage:  python3 extract_legros2019_fig1.py /path/to/1805.02512.pdf [outdir]
Get the PDF:  curl -L -o 1805.02512.pdf https://arxiv.org/pdf/1805.02512v1   (the e-print is the same PDF)
Requires PyMuPDF.

Calibration from the drawn tick marks (PDF points):
  panel c: x = 131.40 (0 K) ... 272.88 (40 K);   y = 387.72 (0) ... 329.64 (15 uOhm cm)
  panel b: x = 363.24 (0 K) ... 512.40 (30 K);   y = 220.44 (0) ... 100.92 (30 uOhm cm)
Drawing coordinates are quantized to 0.12 pt, i.e. 0.034 K and 0.031 uOhm cm in panel c.
"""
import sys, os, json
import numpy as np
import pymupdf


def same(c, ref, tol=0.01):
    return c is not None and all(abs(a - b) < tol for a, b in zip(c, ref))


def main(pdf, outdir):
    page = pymupdf.open(pdf)[10]
    dr = page.get_drawings()
    # panel c: the only blue stroke path with > 1000 segments
    pts = []
    for d in dr:
        if same(d.get('color'), (0, 0, 1)) and len(d['items']) > 1000:
            for it in d['items']:
                if it[0] == 'l':
                    pts += [(it[1].x, it[1].y), (it[2].x, it[2].y)]
    a = np.array(pts)
    T = (a[:, 0] - 131.40) * 40.0 / (272.88 - 131.40)
    r = (387.72 - a[:, 1]) * 15.0 / (387.72 - 329.64)
    o = np.argsort(T, kind='stable'); T, r = T[o], r[o]
    # panel b: orange filled markers (diameter 3.36-3.48 pt) inside the frame; the legend glyphs are excluded by size
    dots = set()
    for d in dr:
        r_ = d['rect']
        if (same(d.get('fill'), (1, 0.55, 0)) and 3.3 < r_.width < 3.5 and 363 < r_.x0 < 522 and r_.y1 < 222):
            dots.add((round((r_.x0 + r_.x1) / 2, 2), round((r_.y0 + r_.y1) / 2, 2)))
    b = np.array(sorted(dots))
    Tb = (b[:, 0] - 363.24) * 30.0 / (512.40 - 363.24)
    rb = (220.44 - b[:, 1]) * 30.0 / (220.44 - 100.92)
    out = dict(
        NdLSCO_p024_33T=dict(T=T.tolist(), drho=r.tolist(),
                             source='Legros et al. 2019 Fig. 1c (data of Daou et al. 2009)', field_T=33),
        LSCO_p023_48T=dict(T=Tb.tolist(), drho=rb.tolist(),
                           source='Legros et al. 2019 Fig. 1b (data of Cooper et al. 2009)', field_T=48),
        units='T in K; drho = rho(T) - rho_0 in micro-Ohm cm')
    os.makedirs(outdir, exist_ok=True)
    json.dump(out, open(os.path.join(outdir, 'legros2019_fig1_rho.json'), 'w'))
    print(f"Nd-LSCO p=0.24, 33 T: {len(T)} points, T = {T.min():.2f}-{T.max():.2f} K, drho up to {r.max():.2f}")
    print(f"LSCO p=0.23, 48 T: {len(Tb)} points, T = {Tb.min():.2f}-{Tb.max():.2f} K")


if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2] if len(sys.argv) > 2 else os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'data'))
