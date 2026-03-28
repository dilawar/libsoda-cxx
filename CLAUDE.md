# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

`libsoda-cxx` is a C++17 port of the LSODA ODE solver library (originally from Heng Li's C implementation). LSODA automatically switches between Adams (non-stiff) and BDF (stiff) methods with adaptive step size and order selection.

## Build Commands

```bash
# Configure and build
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build

# Run tests
cd build && ctest

# Run a single test executable directly
./build/test_lsoda

# Run benchmarks
./build/benchmark_lsoda

# Using just (runs cmake configure + build + ctest)
just
```

## Code Formatting

```bash
clang-format -i src/LSODA.cpp src/LSODA.h
```

Style: WebKit-based, 80-column limit, 4-space indent (see `.clang-format`).

## Architecture

The library consists of two main files:

- **`src/LSODA.h`** — Class definition with the `LSODA` class and `LSODA_ODE_SYSTEM_TYPE` function pointer typedef
- **`src/LSODA.cpp`** — Full solver implementation (~2100 lines)
- **`src/helper.h`** — Template utilities for float comparison (`areEqual`) and vector printing

### Key API

```cpp
// ODE system function signature
typedef void (*LSODA_ODE_SYSTEM_TYPE)(double t, double* y, double* dydt, void* data);

LSODA solver;

// Low-level interface
solver.lsoda(f, neq, y, &t, tout, itol, rtol, atol, itask, &istate, iopt, jt, data);

// Higher-level stepping interface
solver.lsoda_update(f, neq, y, yout, &t, tout, &istate, data, rtol, atol);
```

`istate` starts at 1; successful steps return 2; errors return negative values. Reset to 1 to restart integration.

### Internal Structure

The `LSODA` class holds all integration state as private members — this means **one instance per integration thread** (the benchmark demonstrates thread-safe parallel usage by creating a new `LSODA` instance per thread).

Method selection is tracked by `meth_` (1=Adams, 2=BDF). The solver maintains coefficient tables `elco`/`tesco`/`el` and work arrays `yh_`, `wm_`, `ewt`, `savf`, `acor`.

## Tests

Tests are in `tests/test_LSODA.cpp` and currently cover:
- Chemical kinetics (`fex`) — 3 equations, stiff
- Van der Pol oscillator (`system_scipy`) — 2 equations, mu=1e4
- GitHub issue #10 system (`system_github_issue_10`) — stiff with trigonometric forcing
- Exponential decay (`exponential_decay_rhs`) — 1 equation, analytical solution
- Simple harmonic oscillator (`sho_rhs`) — 2 equations, analytical solution
- Coupled decay chain (`coupled_decay_rhs`) — 3 equations, Bateman equations
- Rational ODE (`rational_ode_rhs`) — 2 equations, analytical solution

Assertions use the `CHECK` macro (throws `runtime_error`) which is active in both Debug and Release builds, unlike `assert()` which is disabled by `-DNDEBUG` in Release.

## CI

GitHub Actions (`.github/workflows/cmake-multi-platform.yml`) runs on push/PR to `master` across Ubuntu (GCC + Clang) and Windows (MSVC), Release build only.

## Changelog

Generated with `git-cliff`. Run `git cliff` to regenerate `CHANGELOG.md` from conventional commits.
