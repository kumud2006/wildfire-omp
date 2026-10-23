#!/usr/bin/env bash
# Shared helpers, sourced by every script in this folder.
#
# Environment variables you can override (examples in README):
#   THREADS="1 2 4 8"   thread counts to test  (default: 1,2,4,... up to nproc)
#   REPEATS=3           repetitions per configuration (median is used in plots)
#   STEPS=200           time steps per run
#   SEED=42             RNG seed
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BIN="${BIN:-$ROOT/wildfire}"
RESULTS="$ROOT/results"
mkdir -p "$RESULTS"

REPEATS="${REPEATS:-3}"
STEPS="${STEPS:-200}"
SEED="${SEED:-42}"

if [ ! -x "$BIN" ]; then echo "Binary $BIN not found. Run 'make' first." >&2; exit 1; fi

# Number of logical CPUs. NOTE: if your CPU has hyper-threading, the number of
# PHYSICAL cores is usually half of this. Set THREADS yourself to be safe.
NCORES="$( (nproc 2>/dev/null || sysctl -n hw.ncpu 2>/dev/null || echo 1) )"

if [ -z "${THREADS:-}" ]; then
  THREADS=""; t=1
  while [ "$t" -lt "$NCORES" ]; do THREADS="$THREADS $t"; t=$((t*2)); done
  THREADS="$THREADS $NCORES"
fi
MAXT="$(echo $THREADS | tr ' ' '\n' | sort -n | tail -1)"

# Thread binding is OFF unless you ask for it (bind_study.sh compares both).
export OMP_PROC_BIND="${OMP_PROC_BIND:-false}"
[ "$OMP_PROC_BIND" != "false" ] && export OMP_PLACES="${OMP_PLACES:-cores}"

ensure_header() {                       # $1 = csv file
  [ -s "$1" ] || "$BIN" --csv-header > "$1"
}

# run_case OUTFILE [wildfire args...]  -> appends REPEATS lines to OUTFILE
run_case() {
  local out="$1"; shift
  ensure_header "$out"
  local r
  for r in $(seq 1 "$REPEATS"); do
    "$BIN" --quiet --steps "$STEPS" --seed "$SEED" "$@" >> "$out"
  done
}

echo "[config] BIN=$BIN  THREADS='$THREADS'  REPEATS=$REPEATS  STEPS=$STEPS  NCORES(logical)=$NCORES  BIND=$OMP_PROC_BIND"
