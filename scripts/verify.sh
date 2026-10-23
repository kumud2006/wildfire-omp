#!/usr/bin/env bash
# CORRECTNESS TEST  (Srikripa)
# Runs every kernel / terrain / thread-count / schedule combination alongside the
# serial reference and compares the grid at EVERY time step (FNV checksum).
# Prints PASS/FAIL per case and exits with an error if anything fails.
source "$(dirname "$0")/common.sh"

SIZES="${SIZES:-100 257 512}"          # odd sizes on purpose (tests boundary handling)
VSTEPS="${VSTEPS:-100}"
fail=0; total=0
printf "%-7s %-6s %-5s %-4s %-8s %-6s %s\n" kernel terrain size thr sched chunk result
for n in $SIZES; do
 for terrain in flat hilly; do
  for kernel in base opt; do
   for th in $THREADS; do
    for cfg in "static 0" "dynamic 4" "guided 8"; do
      set -- $cfg
      total=$((total+1))
      if out=$("$BIN" --quiet --verify --kernel "$kernel" --terrain "$terrain" --n "$n" \
                      --steps "$VSTEPS" --threads "$th" --sched "$1" --chunk "$2" \
                      --ignitions 3 --seed "$SEED" 2>/dev/null); then res=PASS; else res=FAIL; fail=$((fail+1)); fi
      printf "%-7s %-6s %-5s %-4s %-8s %-6s %s\n" "$kernel" "$terrain" "$n" "$th" "$1" "$2" "$res"
    done
   done
  done
 done
done
echo "----------------------------------------------"
echo "$((total-fail)) / $total cases PASSED"
[ "$fail" -eq 0 ] && echo "ALL CORRECT: parallel == serial for every case." || { echo "FAILURES FOUND"; exit 1; }
