# Parallel Wildfire Simulation with OpenMP

A 2D wildfire simulation written in C and parallelized with OpenMP. Fire spreads across a grid based on fuel, wind, and terrain slope. The project compares serial and parallel execution, scheduling choices, SIMD, and scaling.

![Wildfire simulation](results/demo_fire.png)

## Features

- Serial and OpenMP simulation kernels
- Configurable grid size, wind, terrain, and thread count
- Correctness checks against the serial version
- Benchmark scripts and plots for scheduling, SIMD, and scaling

## Run

Requires GCC with OpenMP support. On Linux or WSL:

```bash
make
./wildfire --n 64 --m 64 --steps 20 --kernel opt --threads 4 --verify --ascii
