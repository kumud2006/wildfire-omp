# Build:  make            -> ./wildfire        (optimised, SIMD on)
#         make novec      -> ./wildfire_novec  (SIMD + auto-vectorisation OFF, for comparison)
#         make vecreport  -> shows which loops the compiler vectorised
#         make clean
#
# -ffp-contract=off : stops the compiler from fusing a*b+c differently in
#                     different loops, so serial / base / opt give IDENTICAL
#                     results (needed for the checksum verification).
CC      ?= gcc
CFLAGS  = -std=c11 -O3 -march=native -fopenmp -ffp-contract=off -Wall -Wextra
LDFLAGS = -lm
SRC     = src/main.c src/sim.c src/init.c
HDR     = src/wildfire.h

all: wildfire

wildfire: $(SRC) $(HDR)
	$(CC) $(CFLAGS) $(SRC) -o $@ $(LDFLAGS)

novec: wildfire_novec
wildfire_novec: $(SRC) $(HDR)
	$(CC) $(CFLAGS) -DNO_SIMD -fno-tree-vectorize -fno-tree-slp-vectorize $(SRC) -o $@ $(LDFLAGS)

vecreport:
	$(CC) $(CFLAGS) -fopt-info-vec-optimized -c src/sim.c -o /dev/null 2>&1 | grep -i "sim.c" || true

clean:
	rm -f wildfire wildfire_novec

.PHONY: all novec vecreport clean
