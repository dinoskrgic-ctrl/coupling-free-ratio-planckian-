# Reproducibility package for the v6 Letter
# Claude Opus 5.5 assisted with the code writing and repository formation
This package is for "Planckian dissipation per unit critical entropy as a coupling-free probe of the strange-metal bath":

Manuscript currently in submission, will post arXiv link once it's out 

The script regenerates every computed number and every figure used in the work. Numbers that were read from published papers, rather than computed, are listed under [Literature inputs](#literature-inputs) with their sources.

## Quick start

```bash
pip install numpy scipy matplotlib SciencePlots pymupdf   # pymupdf only for DOWNLOAD=1
cd repro
bash run_all.sh                 # about 50 s; uses the curves already extracted into data/
DOWNLOAD=1 bash run_all.sh      # first fetches the three arXiv files and re-extracts every published curve
FULL=1 bash run_all.sh          # also re-solves the two-band SYK model (several CPU hours, 1-2 GB RAM)
```

Call the scripts through `bash`, not `./`. On some shared or network mounts the execute bit is missing.

Every Python step runs as `python3 -I -B`:
- `-I` means the user site and `PYTHONPATH` are ignored. Each script adds its own folder to `sys.path`.
- `-B` means no `__pycache__` is written.

### What a run writes

- `outputs/*.txt`
- `data/universal_curves.json`
- `data/drift_summary.json`
- `data/michon2019_B_fits.json`
- The figures, as PDF and PNG, in the folder above this one, next to the `.tex` files.

A rerun on 2026-10-21 reproduced `outputs/v6_numbers.txt`, `drift_test.txt`, `band_dos.txt` and `figS_syk_table.txt` identically.

### Tested with

| Package | Version |
|---|---|
| Python | 3.13.16 |
| numpy | 2.5.3 (≥ 2.0 is needed, for `np.trapezoid`) |
| scipy | 1.18.1 |
| matplotlib | 3.11.2 |
| SciencePlots | 2.2.2 |
| PyMuPDF | 2.0 |

### Figure style

The figures use the SciencePlots style `["science", "no-latex"]`, so they need no TeX installation. If you have LaTeX, set `STYLE = ["science"]` at the top of `scripts/make_figures_v6.py`. The fonts will then match the paper exactly.

All figures are single-column (3.375 in wide) and carry no titles; the captions in the `.tex` files describe them.

## Layout

```
repro/
  run_all.sh                     the whole pipeline
  scripts/
    bath.py                      shared bath functionals: beta = S(0), C_b, arctan table, optical model
    extract_michon2023_fig2.py   Michon et al. 2023 Figs. 2a (rho) and 2b (optical rate) from the arXiv source fig02.pdf
    extract_legros2019_fig1.py   Legros et al. 2019 Fig. 1b (LSCO p=0.23, 48 T) and 1c (Nd-LSCO p=0.24, 33 T)
    extract_michon2019_curves.py Michon et al. 2019 Figs. 2d, 3d (C_el/T of Eu- and Nd-LSCO)
    fit_michon2019_B.py          logarithmic fits C/T = B ln(T0/T)
    universal_thermal_mass.py    large-M thermal-mass saddle of the self-damped field (beta, C_b, delta, drift)
    band_dos.py                  band systematic of gamma_b and (m/m_e)K for LSCO-type bands at n = 0.76
    drift_test_v6.py             predicted resistivity exponent n vs power-law fits to three published curves
    v6_numbers.py                every computed number, by section -> outputs/v6_numbers.txt
    make_figures_v6.py           fig1_map, fig2_field, fig3_data, figS_spectra, figS_exponents, figS_materials, figS_drift
    make_figS_syk_v6.py          figS_syk and outputs/figS_syk_table.txt
  syk_two_band/                  large-N two-band SYK solver (Matsubara free energy + real-axis Kubo rate)
    run_syk.sh                   --quick (default) re-fits the saved out/*.npy arrays; --full re-solves (hours)
    out/                         saved free energies F_*.npy, Kubo rates real_alpha*.npy, solver logs
  data/                          extracted curves and computed tables (provenance below)
  outputs/                       text outputs of the last run
  v4_checks/                     (zip only) the archived v4 package, for the inherited numbers listed below
```

## Where each number comes from

### Letter

| Item | Script | Output |
|---|---|---|
| C_b ≃ 3.6β (3.48–3.76), R_opt = 2C_b/π | `v6_numbers.py` (uses `bath.py`) | `v6_numbers.txt` § S2–S3 |
| Q_conf = 4.913, b = 3/2 | `v6_numbers.py` | § S4 |
| Two-band SYK: Q/Q_conf = 1.0020–1.0024 | `syk_two_band/` + `make_figS_syk_v6.py` | `figS_syk_table.txt` |
| Field: δ = 0.007–0.074, η = 0.93–0.99; floor β ≥ 2.8 (2.8–3.8); Q ≈ 9–13 (up to 15.5); drift 0.1–0.5 per decade | `universal_thermal_mass.py` | `universal_thermal_mass.txt`, `v6_numbers.txt` § S5 |
| Field identity B = γ_b λ | algebra, SM Eq. (S8) | `v6_numbers.txt` § S5 (stated) |
| LSCO slope 0.616–0.635, λ = 0.40 | `extract_michon2023_fig2.py`, `v6_numbers.py` | § S6 |
| LSCO R_opt 8.0–9.1; C_b, r and β by both routes | `v6_numbers.py` | § S6 |
| LSCO Q = 8.8–10.1 (7.8–11.7); η = 0.70 / 0.81 | `v6_numbers.py` | § S6 |
| Band factors 0.87–2.0 | `band_dos.py` | `band_dos.txt`, `v6_numbers.txt` § S6 |
| Eu-LSCO B = 2.63, power law T^-0.17; Nd-LSCO B = 2.87 above 2 K | `fit_michon2019_B.py` | `fit_michon2019_B.txt`, `v6_numbers.txt` § S7 |
| Nd-LSCO η, both routes and conformal | `v6_numbers.py` | § S7 |
| Drift: predicted n and the three fits | `drift_test_v6.py` | `drift_test.txt`, `v6_numbers.txt` § S8 |
| Figs. 1–3 | `make_figures_v6.py` | `../fig1_map`, `../fig2_field`, `../fig3_data` |

### Supplemental Material

Most SM numbers come from the same scripts and outputs as the Letter. The following come from the archived v4 package, which is unchanged and is included in the zip as `v4_checks/`:

| SM item | Archived script (in `v4_checks/`) |
|---|---|
| S4: digamma-kernel log slopes 9.873 and 6.579 | `conformal/entropy_dissipation_ratio.py` |
| S5: reservoir route, β = 3.35 / 3.26 / 3.00 / 2.69 at η = 1 / 0.5 / 0.1 / 0.02 | `toy_models/toy_two_damping.py` |
| S7: Nd-LSCO Q at B = 3.23 (7.41 ± 1.22 in-plane, 4.07 ± 1.45 ADMR), rescaled by `v6_numbers.py` to B = 2.87 | `conformal/stress_tests.py`, `conformal/entropy_dissipation_ratio.py` |
| S7: Sr₃Ru₂O₇ Q = 4.7 (3.7–6.1) | `q_survey/sr327_band.py` |
| S7: Kondo-lattice Q values | `q_survey/q_survey.py` |
| S11: Lorenz ratios 0.713, 0.76–0.80, 0.731 and 0.82 | `shapes/lorenz_shape.py`, `added_checks/manuscript_numbers.py` |
| S11: χ''/ω half point 2.33 r | `added_checks/v4_numbers.py`, item (e) |
| S12: spin fluctuations at r = 2.9 give R_opt ≈ 2.3 | `added_checks/v4_numbers.py`, item (f) |
| S12: pair-fluctuation share ≤ 30% | `failed_attempts/pair_fluctuations.py` |

Run them with `bash v4_checks/run_all.sh`. Its SYK step reports FAIL because `syk_two_band/` lives at the top level of this package; the rest runs.

That archive's `data/michon2023_fig2b_curves.json` is the old extraction. Its zero sits 0.91 too high in ħ/τk_BT (see below). Only slopes were ever taken from it. The corrected file is `data/michon2023_fig2b_curves.json` in this package.

### Not scripted

The CeCoIn₅ estimate b ≈ 17 (13–29), or 2–7 under other input choices (SM S7), is a hand estimate. It uses the terahertz data of Pracht et al. (arXiv:2604.10621), B at 5 T and γ_b = 3.3 mJ K⁻² mol⁻¹. The working is in the project note `bath-realizations/12_two_b_tests.md`.

## Data provenance

| File | Source | Extraction and checks |
|---|---|---|
| `data/michon2023_fig2b_curves.json` | Michon et al., *Nat. Commun.* **14**, 3033 (2023), arXiv:2205.04030, Fig. 2b | See [the Michon 2023 extraction](#the-michon-2023-extraction) below |
| `data/michon2023_fig2a_rho.json` | same paper, Fig. 2a (ρ at 0 and 16 T) | Same calibration as Fig. 2b. Check: the dashed guide line comes out as ρ = 13.5 + 0.629 T µΩ cm, against the quoted slope 0.63 |
| `data/legros2019_fig1_rho.json` | Legros et al., *Nat. Phys.* **15**, 142 (2019), arXiv:1805.02512v1, page index 10 | See [the Legros 2019 extraction](#the-legros-2019-extraction) below |
| `data/michon2019_curves.json` | Michon et al., *Nature* **567**, 218 (2019), arXiv:1804.08502, Figs. 2d and 3d (page indices 9 and 11) | Stroke points selected by colour, with axes calibrated at y = 0 and 25 mJ K⁻² mol⁻¹ and x = 1 and 10 K. Check: the reference lines come out at 4.14 and 3.62, against 4.2 and 3.6 in the paper. About 3,700 points per curve |
| `data/michon2019_B_fits.json` | derived | `fit_michon2019_B.py` |
| `data/universal_curves.json` | computed | `universal_thermal_mass.py` (no experimental input) |
| `data/drift_summary.json` | derived | `drift_test_v6.py` |
| `syk_two_band/out/*` | computed | SYK solvers. The 2²¹-point real-axis self-energies `Sc_*.npy` are not included: they are written by `--full` but not needed to re-fit |

### The Michon 2023 extraction

- The curves are vector paths of the arXiv source file `fig02.pdf`: 8 curves, 50–300 K, identified by legend colour. The 40 K curve leaves the frame.
- The axes are calibrated on the frame lines and major tick marks: x from 263.97 to 434.04 pt for 0–100, and y from 176.80 to 6.73 pt for 0–70.
- **An earlier extraction, whose code was lost, put the zero at the centre of the "0" label.** That is 2.2 pt below the axis, which made every rate 0.91 too high. Slopes were unaffected.

### The Legros 2019 extraction

- **Panel c** (Nd-LSCO, p = 0.24, 33 T; Daou et al. 2009) is a vector curve of 3,246 points. It is calibrated on tick marks: x from 131.40 to 272.88 pt for 0–40 K, and y from 387.72 to 329.64 pt for 0–15 µΩ cm.
- **Panel b** (LSCO, p = 0.23, 48 T; Cooper et al. 2009) is 16 orange dots of width 3.3–3.5 pt. It is calibrated with x from 363.24 to 512.40 pt for 0–30 K, and y from 220.44 to 100.92 pt for 0–30 µΩ cm.

## Literature inputs

These are read from the text or figures of published papers and typed into the scripts. They are in `INPUTS` in `scripts/v6_numbers.py` and at the top of `scripts/make_figures_v6.py`.

| Input | Value | Source |
|---|---|---|
| Band Sommerfeld coefficient of LSCO at p = 0.24 | γ_b = 3.89 mJ K⁻² mol⁻¹ (t′/t = −0.17, m_b = 2.76 m_e) | Michon 2023, tight-binding band |
| Zero-frequency rate, optics | ħ/τ(0)k_BT = 5 | Michon 2023, text |
| Zero-frequency rate, dc | 5.7 | Michon 2023: A = 0.63 µΩ cm/K with their conversion 5 ↔ 0.55 |
| Coupling from the optical mass logarithm | λ = 0.46 | 2 × Michon 2023's g = 0.23 |
| Logarithmic coefficient of LSCO | B = 2.2 ± 0.3 mJ K⁻² mol⁻¹ | Girod et al. 2021 |
| Release ratio | 1.47–1.66 | α / (C_b,conf × 0.46) with α = 5–5.7; the same as C_b = α/λ at λ = 0.46 |
| Nd-LSCO in-plane slope | A₁□ = 7.4 ± 0.8 Ω/K | Daou 2009 / Legros 2019 |
| Nd-LSCO ADMR rate and band | ADMR rate and band as published | Grissonnanche et al. 2021 (via the archived `stress_tests.py`) |
| Nodal photoemission exponent | β = 3–4 | Reber et al. 2019; Smit et al. 2021 (Bi2212, Bi2201) |
| Other materials | Sr₃Ru₂O₇ and Kondo-lattice inputs | Sources listed in `v4_checks/README.md` ("Experimental inputs typed directly into scripts") |

## Known limits of this package

- **The SYK solvers were not re-run here.** The quick path re-fits the saved free energies and Kubo rates; `FULL=1` re-solves them.
- **The extraction scripts depend on the arXiv file versions named above.** They use hard-coded page indices and point coordinates. A different arXiv version may need new calibration numbers, which the scripts print checks for.

## Changes for the combined single-file version (2026-10-22)

- `drift_test_v6.py` also fits the zero-field LSCO curve over the predicted window, 40–85 K, and its sub-windows
  (n = 1.031, 1.022–1.039), and prints how far the 48 T LSCO fit sits above the predicted range (1.4 and 2.3 s.d.).
  These are the numbers in the Letter's drift paragraph and in Sec. S8.
- `make_figures_v6.py` draws the zero-field LSCO box in Fig. 2b over 40–85 K.
- `make_fig_lattice.py` draws the first-page lattice schematic, `fig_lattice.pdf`.
