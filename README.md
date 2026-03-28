[![CMake on multiple platforms](https://github.com/dilawar/libsoda-cxx/actions/workflows/cmake-multi-platform.yml/badge.svg)](https://github.com/dilawar/libsoda-cxx/actions/workflows/cmake-multi-platform.yml)

# libsoda++

C++17 port of the LSODA ODE solver. LSODA automatically switches between
Adams (non-stiff) and BDF (stiff) integration methods with adaptive step size.

## Usage

```cpp
#include "LSODA.h"

void my_ode(double t, double* y, double* dydt, void* data) {
    dydt[0] = -y[0];  // dy/dt = -y
}

LSODA solver;
std::vector<double> y = {1.0}, yout;
double t = 0.0;
int istate = 1;
solver.lsoda_update(my_ode, 1, y, yout, &t, 1.0, &istate, nullptr);
// yout[1] holds the result (1-indexed)
```

See [tests/test_LSODA.cpp](tests/test_LSODA.cpp) for more examples.

## Build

```sh
cmake -S . -B build && cmake --build build && cd build && ctest
```

## Credits

Based on Heng Li's C port: <http://lh3lh3.users.sourceforge.net/download/lsoda.c>

## Related

- https://github.com/sdwfrost/liblsoda
