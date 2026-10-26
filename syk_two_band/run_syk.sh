#!/usr/bin/env bash
# Two-band SYK check (v6 Supplemental Material Sec. S10, Q/Q_conf = 1.0020-1.0024).
#   bash run_syk.sh --quick   (default) re-fit the saved out/*.npy arrays: analyze.py + ../scripts/make_figS_syk_v6.py (seconds)
#   bash run_syk.sh --full    re-solve everything (Matsubara free energies + real-axis Kubo rates on 2^21-point grids),
#                          then the same two steps. Expect several hours on one core and ~1-2 GB RAM; the original
#                          real-axis logs show 90-270 s per temperature point, 10 points per run, 3 real-axis runs.
# The full run OVERWRITES out/*.npy and out/log_*.txt in this folder.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PY="python3 -I -B"   # -I: isolated (no user site, script dir not on sys.path); -B: no __pycache__
MODE="${1:---quick}"
if [[ "$MODE" == "--full" ]]; then
  O="$HERE/out"
  # Matsubara free energies (frozen bath, Uf = 0.05): F_U{U}_dt0.02.npy
  for U in 0.17 0.3 0.5; do $PY "$HERE/run_matsu.py" $U 0.02 > "$O/log_U${U}_0.02.txt"; done
  # variants: half time step; Uf = 0.1; bath feedback N_c = N_f
  $PY "$HERE/run_matsu2.py" 0.3 0.01 0.0 0.05 > "$O/log_dt01.txt"
  $PY "$HERE/run_matsu2.py" 0.3 0.02 0.0 0.1  > "$O/log_uf01.txt"
  $PY "$HERE/run_matsu2.py" 0.3 0.02 1.0 0.05 > "$O/log_r1.txt"
  # real-axis Kubo rates
  $PY "$HERE/run_real.py"        > "$O/log_real.txt"
  $PY "$HERE/run_more.py" uf0.1  > "$O/log_real_uf01.txt"
  $PY "$HERE/run_more.py" fc     > "$O/log_real_fc.txt"
elif [[ "$MODE" != "--quick" ]]; then
  echo "usage: $0 [--quick|--full]" >&2; exit 2
fi
echo "== analyze.py (B fits of the saved free energies)"
$PY "$HERE/analyze.py"
echo "== make_figS_syk_v6.py (figS_syk.{pdf,png} next to the .tex files; ../outputs/figS_syk_table.txt)"
$PY "$HERE/../scripts/make_figS_syk_v6.py"
