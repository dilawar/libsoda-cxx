#include <lsoda/LSODA.h>
#include <cmath>
#include <stdexcept>
#include <vector>

// Simple exponential decay: dy/dt = -y, y(0)=1, solution: y(t)=exp(-t)
static void exponential_decay(double, double* y, double* dydt, void*)
{
    dydt[0] = -y[0];
}

int main()
{
    LSODA solver;
    std::vector<double> y = { 1.0 }, yout;
    double t = 0.0;
    int istate = 1;

    solver.lsoda_update(exponential_decay, 1, y, yout, &t, 1.0, &istate, nullptr, 1e-8, 1e-10);

    if (istate <= 0)
        throw std::runtime_error("lsoda_update failed: istate=" + std::to_string(istate));

    double expected = std::exp(-1.0);
    double err = std::fabs(yout[1] - expected);
    if (err > 1e-6)
        throw std::runtime_error(
            "result out of tolerance: got " + std::to_string(yout[1])
            + ", expected " + std::to_string(expected));

    return 0;
}
