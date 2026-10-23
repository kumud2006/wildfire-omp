#!/usr/bin/env bash
# THREAD BINDING STUDY  (Kumud): OMP_PROC_BIND=false vs close vs spread.
# Output: results/bind.csv   (the 'tag' column holds the binding mode)
source "$(dirname "$0")/common.sh"
OUT="$RESULTS/bind.csv"
N="${N:-2048}"
TH="${TH:-$MAXT}"
for mode in false close spread; do
  echo "== OMP_PROC_BIND=$mode"
  if [ "$mode" = false ]; then unset OMP_PLACES; else export OMP_PLACES=cores; fi
  OMP_PROC_BIND=$mode run_case "$OUT" --kernel opt --terrain hilly --n "$N" --threads "$TH" --tag "bind-$mode"
done
echo "Done -> $OUT"
