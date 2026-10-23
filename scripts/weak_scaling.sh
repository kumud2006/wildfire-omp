#!/usr/bin/env bash
# WEAK SCALING  (Gautami): work per thread stays constant -> grid area grows with
# the thread count:  N(p) = BASE_N * sqrt(p).   Ideal: run time stays flat.
# Output: results/weak.csv
source "$(dirname "$0")/common.sh"
OUT="$RESULTS/weak.csv"
BASE_N="${BASE_N:-1024}"
KERNELS="${KERNELS:-opt}"
TERRAIN="${TERRAIN:-hilly}"

for k in $KERNELS; do
  for th in $THREADS; do
    n=$(awk -v b="$BASE_N" -v p="$th" 'BEGIN{printf "%d", b*sqrt(p)+0.5}')
    echo "== kernel=$k threads=$th N=$n"
    run_case "$OUT" --kernel "$k" --terrain "$TERRAIN" --n "$n" --threads "$th"
  done
done
echo "Done -> $OUT"
