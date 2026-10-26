"""Lattice schematic for the first page (from the left panel of the v4 Fig. 1, same drawing and random seed).

A two-dimensional lattice of nearly fixed-amplitude local units (arrows), critical at p*, sitting in the
conduction-electron sea (grey wavy lines) and coupled to the electron density only through spatially random
local fields of zero mean (red +, blue -). No text in the figure; everything is explained in the caption.

Run: python3 make_fig_lattice.py   (numpy, matplotlib). Writes fig_lattice.pdf and .png next to the .tex files (../../).
"""
import os
import numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle

HERE = os.path.dirname(os.path.abspath(__file__))
W1 = 3.375                                   # PRL single-column width, inches
POS, NEG = "#d6604d", "#4393c3"              # sign of the random local field


def fig_lattice():
    fig = plt.figure(figsize=(W1, W1*5.3/6.9))
    a = fig.add_axes([0, 0, 1, 1])
    a.set_xlim(-0.75, 6.15); a.set_ylim(-0.45, 4.85); a.set_aspect("equal"); a.axis("off")
    rng = np.random.default_rng(7)
    xs = np.linspace(-0.8, 6.2, 300)
    for j in np.arange(-0.25, 4.6, 0.62):                                   # electron sea
        a.plot(xs, j + 0.06*np.sin(5.5*xs + 3*j), color="0.80", lw=0.7, zorder=0)
    for (x, y) in rng.uniform([-0.4, -0.3], [5.4, 4.3], (11, 2)):          # random fields, zero mean
        s = rng.choice([-1, 1]); c = POS if s > 0 else NEG
        a.add_patch(Circle((x, y), 0.48, fc=c, ec="none", alpha=0.20, zorder=1))
        a.plot(x, y, marker="+" if s > 0 else "_", ms=6, mew=1.3, color=c, zorder=2)
    for i in range(6):                                                      # fixed-amplitude local units
        for j in range(5):
            th = rng.normal(0, 0.35) + (0 if rng.random() < 0.5 else np.pi)
            dx, dy = 0.28*np.sin(th), 0.28*np.cos(th)
            a.annotate("", xy=(i+dx, j+dy), xytext=(i-dx, j-dy), zorder=4,
                       arrowprops=dict(arrowstyle="-|>", color="k", lw=0.9, mutation_scale=6, shrinkA=0, shrinkB=0))
            a.add_patch(Circle((i, j), 0.06, fc="k", ec="none", zorder=5))
    for i in range(6): a.plot([i, i], [0, 4], color="0.6", lw=0.5, zorder=3)
    for j in range(5): a.plot([0, 5], [j, j], color="0.6", lw=0.5, zorder=3)
    fig.savefig(os.path.join(HERE, "..", "..", "fig_lattice.pdf"))
    fig.savefig(os.path.join(HERE, "..", "..", "fig_lattice.png"), dpi=300)
    plt.close(fig)


if __name__ == "__main__":
    fig_lattice()
