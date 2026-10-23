#!/usr/bin/env bash
# SCHEDULING STUDY  (Kumud): static / dynamic / guided, several chunk sizes,
# flat terrain (even work) vs hilly terrain (uneven work).
# Output: results/schedule.csv   (the 'imbalance' column = max/mean thread busy time)
#
# Chunk meaning:  base kernel -> CELLS (collapse(2));  opt kernel -> ROWS.
source "$(dirname "$0")/common.sh"
OUT="$RESULTS/schedule.csv"
N="${N:-2048}"
TH="${TH:-$MAXT}"                                    # thread count for the study
CHUNKS_BASE="${CHUNKS_BASE:-0 64 1024 16384}"        # cells
CHUNKS_OPT="${CHUNKS_OPT:-0 1 4 16}"                 # rows
TERRAINS="${TERRAINS:-flat hilly}"
KERNELS="${KERNELS:-base opt}"

for terrain in $TERRAINS; do
 for k in $KERNELS; do
  if [ "$k" = base ]; then CH="$CHUNKS_BASE"; else CH="$CHUNKS_OPT"; fi
  for sched in static dynamic guided; do
   for c in $CH; do
     echo "== $terrain kernel=$k sched=$sched chunk=$c threads=$TH"
     run_case "$OUT" --kernel "$k" --terrain "$terrain" --n "$N" --threads "$TH" \
                     --sched "$sched" --chunk "$c" --ignitions 4
   done
  done
 done
done
echo "Done -> $OUT"
