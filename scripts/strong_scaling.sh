#!/usr/bin/env bash
# STRONG SCALING  (Srikripa): fixed problem size, more threads.
# Output: results/strong.csv   (serial baseline rows included)
source "$(dirname "$0")/common.sh"
OUT="$RESULTS/strong.csv"
SIZES="${SIZES:-1024 2048 4096}"
KERNELS="${KERNELS:-base opt}"
TERRAIN="${TERRAIN:-hilly}"

for n in $SIZES; do
  echo "== N=$n : serial baseline"
  run_case "$OUT" --kernel serial --terrain "$TERRAIN" --n "$n"
  for k in $KERNELS; do
    for th in $THREADS; do
      echo "== N=$n kernel=$k threads=$th"
      run_case "$OUT" --kernel "$k" --terrain "$TERRAIN" --n "$n" --threads "$th"
    done
  done
done
echo "Done -> $OUT"
