#!/usr/bin/env bash
# Regenerate every computed number and every figure of the v6 Letter (entropy_prl_v6.tex) and its
# Supplemental Material (entropy_prl_v6_SM.tex). See README.md for what each step does and where its inputs come from.
#
#   bash run_all.sh                 default: uses the curves already extracted into data/ (about a minute)
#   DOWNLOAD=1 bash run_all.sh      first fetches the three arXiv files and re-extracts the published curves (needs
#                                   network access and PyMuPDF); the downloads go to downloads/, one folder per paper
#   FULL=1 bash run_all.sh          also re-solves the two-band SYK model from scratch (several CPU hours)
#
# Outputs: outputs/*.txt, data/universal_curves.json, data/drift_summary.json, and the figures (PDF + PNG) in the
# folder above this one, next to the .tex files. Call through bash (the shared mount may not allow ./run_all.sh).
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
S="$HERE/scripts"
PY="python3 -I -B"   # -I: no user site, script dir not on sys.path (scripts add their own folder); -B: no __pycache__
mkdir -p "$HERE/outputs"

step() { echo; echo "=== $1"; }

if [[ "${DOWNLOAD:-0}" == "1" ]]; then
  DL="$HERE/downloads"
  mkdir -p "$DL/michon2023" "$DL/legros2019" "$DL/michon2019"
  step "download arXiv:2205.04030 source (Michon et al. 2023; Fig. 2 is fig02.pdf)"
  curl -sSL -o "$DL/michon2023/src.tar" https://arxiv.org/e-print/2205.04030
  tar -xf "$DL/michon2023/src.tar" -C "$DL/michon2023" fig02.pdf
  step "download arXiv:1805.02512v1 (Legros et al. 2019; Fig. 1 on page index 10)"
  curl -sSL -o "$DL/legros2019/1805.02512v1.pdf" https://arxiv.org/pdf/1805.02512v1
  step "download arXiv:1804.08502 (Michon et al. 2019; Figs. 2d and 3d on page indices 9 and 11)"
  curl -sSL -o "$DL/michon2019/1804.08502.pdf" https://arxiv.org/pdf/1804.08502
  step "extract Michon 2023 Figs. 2a and 2b -> data/michon2023_fig2{a_rho,b_curves}.json"
  $PY "$S/extract_michon2023_fig2.py" "$DL/michon2023/fig02.pdf" "$HERE/data"
  step "extract Legros 2019 Fig. 1b,c -> data/legros2019_fig1_rho.json"
  $PY "$S/extract_legros2019_fig1.py" "$DL/legros2019/1805.02512v1.pdf" "$HERE/data"
  step "extract Michon 2019 C/T curves -> data/michon2019_curves.json"
  $PY "$S/extract_michon2019_curves.py" "$DL/michon2019/1804.08502.pdf" "$HERE/data/michon2019_curves.json"
  FIT_PDF="$DL/michon2019/1804.08502.pdf"
else
  FIT_PDF="-"
fi

step "fit the C/T logarithms (B) -> data/michon2019_B_fits.json"
$PY "$S/fit_michon2019_B.py" "$FIT_PDF" "$HERE/data/michon2019_curves.json" "$HERE/data/michon2019_B_fits.json" \
  | tee "$HERE/outputs/fit_michon2019_B.txt"
step "thermal-mass saddle of the self-damped field -> data/universal_curves.json"
$PY "$S/universal_thermal_mass.py" | tee "$HERE/outputs/universal_thermal_mass.txt"
step "band systematic -> outputs/band_dos.txt"
$PY "$S/band_dos.py" | tee "$HERE/outputs/band_dos.txt"
step "drift test -> outputs/drift_test.txt, data/drift_summary.json"
$PY "$S/drift_test_v6.py" > /dev/null
step "number list -> outputs/v6_numbers.txt"
$PY "$S/v6_numbers.py" > /dev/null
step "figures fig1_map, fig2_field, fig3_data, figS_spectra, figS_exponents, figS_materials, figS_drift"
$PY "$S/make_figures_v6.py"
$PY "$S/make_fig_lattice.py"
step "two-band SYK check and figS_syk"
if [[ "${FULL:-0}" == "1" ]]; then
  bash "$HERE/syk_two_band/run_syk.sh" --full | tee "$HERE/outputs/syk_two_band.txt"
else
  bash "$HERE/syk_two_band/run_syk.sh" --quick | tee "$HERE/outputs/syk_two_band.txt"
fi
echo; echo "done: see outputs/v6_numbers.txt"
