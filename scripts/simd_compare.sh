#!/usr/bin/env bash
# SIMD / MEMORY STUDY  (Gautami): compare
#     base  (collapse(2), no tricks)
#     opt   built WITHOUT vectorisation   (./wildfire_novec)
#     opt   built WITH SIMD               (./wildfire)
# Output: results/simd.csv  ('tag' column: base | opt-novec | opt-simd)
source "$(dirname "$0")/common.sh"
OUT="$RESULTS/simd.csv"
NOVEC="$ROOT/wildfire_novec"
[ -x "$NOVEC" ] || { echo "Run 'make novec' first." >&2; exit 1; }
SIZES="${SIZES:-1024 2048}"
TH_LIST="${TH_LIST:-1 $MAXT}"
for n in $SIZES; do
 for th in $TH_LIST; do
   echo "== N=$n threads=$th"
   run_case "$OUT" --kernel base --n "$n" --threads "$th" --tag base
   BIN="$NOVEC" run_case "$OUT" --kernel opt --n "$n" --threads "$th" --tag opt-novec
   run_case "$OUT" --kernel opt --n "$n" --threads "$th" --tag opt-simd
 done
done
echo "Done -> $OUT"
