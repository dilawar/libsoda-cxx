#include <iostream>
#include <iomanip>
#include <vector>
#include <cassert>
#include <chrono>
#include <cmath>
#include <stdexcept>

#include "../src/LSODA.h"
#include "../src/helper.h"

using namespace std;

// Describe the system.
static void fex(double t, double* y, double* ydot, void* data)
{
    (void)t;
    (void)data;

    ydot[0] = 1.0E4 * y[1] * y[2] - .04E0 * y[0];
    // Don't swap ydot[1] and ydot[2]. The order will change and test will fail.
    ydot[2] = 3.0E7 * y[1] * y[1];
    ydot[1] = -1.0 * (ydot[0] + ydot[2]);
}

static void system_scipy(double t, double* y, double* ydot, void* data)
{
    (void)t;
    (void)data;

    double mu = 1E4;
    ydot[0] = y[1];
    ydot[1] = mu * (1 - y[0] * y[0]) * y[1] - y[0];
}

// Simple harmonic oscillator: dy0/dt = y1, dy1/dt = -omega^2 * y0
// omega = 2 => analytical solution: y0(t) = cos(2t) + 0.05*sin(2t)
static void sho_rhs(double t, double* y, double* ydot, void* data)
{
    (void)t;
    (void)data;
    ydot[0] = y[1];
    ydot[1] = -4.0 * y[0];
}

// Radioactive decay chain: 3 species with decay constants lambda = {0.17, 0.23, 0.29}
// dz0/dt = -0.17*z0
// dz1/dt = -0.23*z1 + 0.17*z0
// dz2/dt = -0.29*z2 + 0.23*z1
static void coupled_decay_rhs(double t, double* y, double* ydot, void* data)
{
    (void)t;
    (void)data;
    ydot[0] = -0.17 * y[0];
    ydot[1] = -0.23 * y[1] + 0.17 * y[0];
    ydot[2] = -0.29 * y[2] + 0.23 * y[1];
}

// Exponential decay: dy/dt = -y, solution: y(t) = exp(-t)
static void exponential_decay_rhs(double t, double* y, double* ydot, void* data)
{
    (void)t;
    (void)data;
    ydot[0] = -y[0];
}

// Rational ODE from SciPy test suite (sol_rational):
// dy0/dt = y1 / t
// dy1/dt = y1 * (y0 + 2*y1 - 1) / (t * (y0 - 1))
// Analytical solution: y = (t/(t+10), 10*t/(t+10)^2)
static void rational_ode_rhs(double t, double* y, double* ydot, void* data)
{
    (void)data;
    ydot[0] = y[1] / t;
    ydot[1] = y[1] * (y[0] + 2.0 * y[1] - 1.0) / (t * (y[0] - 1.0));
}

// This system is described here
// https://github.com/sdwfrost/liblsoda/issues/10
static void system_github_issue_10(
    double t, double* y, double* ydot, void* data)
{
    (void)data;

    ydot[0] = 9 * y[0] + 24 * y[1] + 5 * cos(t) - (1 / 3) * sin(t);
    ydot[1] = -24 * y[0] - 51 * y[1] - 95 * cos(t) + (1 / 3) * sin(t);
}

// Verify exponential decay: dy/dt = -y, y(0)=1
// Integrates from t=0 to t=1, expected y(1) = exp(-1)
int test_exponential_decay()
{
    double t = 0.0, tout = 1.0;
    vector<double> y = { 1.0 };
    int istate = 1;
    LSODA lsoda;
    vector<double> yout;

    lsoda.lsoda_update(
        exponential_decay_rhs, 1, y, yout, &t, tout, &istate, nullptr, 1e-8, 1e-10);

    if (istate <= 0) {
        cerr << "error istate = " << istate << endl;
        throw runtime_error("test_exponential_decay failed");
    }

    double expected = exp(-1.0);
    assert(areEqual(expected, yout[1], 1e-6));
    return 0;
}

// Verify simple harmonic oscillator (from SciPy test suite):
// dy0/dt = y1, dy1/dt = -4*y0
// IC: y = {1.0, 0.1}, integrate from t=0 to t=1.09
// Analytical: y0(t) = cos(2t) + 0.05*sin(2t)
//             y1(t) = -2*sin(2t) + 0.1*cos(2t)
int test_simple_harmonic_oscillator()
{
    double t = 0.0, tout = 1.09;
    vector<double> y = { 1.0, 0.1 };
    int istate = 1;
    LSODA lsoda;
    vector<double> yout;

    lsoda.lsoda_update(
        sho_rhs, 2, y, yout, &t, tout, &istate, nullptr, 1e-8, 1e-10);

    if (istate <= 0) {
        cerr << "error istate = " << istate << endl;
        throw runtime_error("test_simple_harmonic_oscillator failed");
    }

    double omega = 2.0;
    double expected_y0 = cos(omega * tout) + 0.1 / omega * sin(omega * tout);
    double expected_y1 = -omega * sin(omega * tout) + 0.1 * cos(omega * tout);

    assert(areEqual(expected_y0, yout[1], 1e-5));
    assert(areEqual(expected_y1, yout[2], 1e-5));
    return 0;
}

// Verify radioactive decay chain (from SciPy CoupledDecay test):
// dz0/dt = -0.17*z0
// dz1/dt = -0.23*z1 + 0.17*z0
// dz2/dt = -0.29*z2 + 0.23*z1
// IC: {5, 7, 13}, integrate from t=0 to t=0.5
// Expected values from Bateman equations (computed analytically)
int test_coupled_decay()
{
    double t = 0.0, tout = 0.5;
    vector<double> y = { 5.0, 7.0, 13.0 };
    int istate = 1;
    LSODA lsoda;
    vector<double> yout;

    lsoda.lsoda_update(
        coupled_decay_rhs, 3, y, yout, &t, tout, &istate, nullptr, 1e-8, 1e-10);

    if (istate <= 0) {
        cerr << "error istate = " << istate << endl;
        throw runtime_error("test_coupled_decay failed");
    }

    // Bateman equations (analytical solution)
    const double lmbd0 = 0.17, lmbd1 = 0.23, lmbd2 = 0.29;
    const double z0_0 = 5.0, z0_1 = 7.0, z0_2 = 13.0;
    const double d10 = lmbd1 - lmbd0; // 0.06
    const double d21 = lmbd2 - lmbd1; // 0.06
    const double d20 = lmbd2 - lmbd0; // 0.12
    const double e0 = exp(-lmbd0 * tout);
    const double e1 = exp(-lmbd1 * tout);
    const double e2 = exp(-lmbd2 * tout);

    double expected_z0 = z0_0 * e0;
    double expected_z1 = z0_1 * e1 + z0_0 * lmbd0 / d10 * (e0 - e1);
    double expected_z2 = z0_2 * e2 + z0_1 * lmbd1 / d21 * (e1 - e2)
        + lmbd1 * lmbd0 * z0_0 / d10
            * (1.0 / d20 * (e0 - e2) - 1.0 / d21 * (e1 - e2));

    assert(areEqual(expected_z0, yout[1], 1e-6));
    assert(areEqual(expected_z1, yout[2], 1e-6));
    assert(areEqual(expected_z2, yout[3], 1e-6));
    return 0;
}

// Verify rational ODE (from SciPy sol_rational test):
// dy0/dt = y1/t, dy1/dt = y1*(y0 + 2*y1 - 1) / (t*(y0 - 1))
// Analytical solution: y = (t/(t+10), 10*t/(t+10)^2)
// Integrates from t=5 (IC = sol_rational(5)) to t=9
int test_rational_ode()
{
    double t = 5.0, tout = 9.0;
    // IC = analytical solution evaluated at t=5
    vector<double> y = { 5.0 / 15.0, 10.0 * 5.0 / (15.0 * 15.0) };
    int istate = 1;
    LSODA lsoda;
    vector<double> yout;

    lsoda.lsoda_update(
        rational_ode_rhs, 2, y, yout, &t, tout, &istate, nullptr, 1e-8, 1e-10);

    if (istate <= 0) {
        cerr << "error istate = " << istate << endl;
        throw runtime_error("test_rational_ode failed");
    }

    double expected_y0 = tout / (tout + 10.0);
    double expected_y1 = 10.0 * tout / ((tout + 10.0) * (tout + 10.0));

    assert(areEqual(expected_y0, yout[1], 1e-5));
    assert(areEqual(expected_y1, yout[2], 1e-5));
    return 0;
}

int test_github_system(void)
{
    // cout << "Running test given
    // https://github.com/sdwfrost/liblsoda/issues/10" << endl;
    double t = 0e0, tout = 0.5;

    vector<double> y = { 4.0 / 3.0, 2.0 / 3.0 };
    int istate = 1;

    LSODA lsoda;

    vector<double> yout;
    vector<double> res;

    for (size_t i = 0; i < 10; i++) {
        lsoda.lsoda_update(
            system_github_issue_10, 2, y, yout, &t, tout, &istate, nullptr);
        res.push_back(yout[1]);
        res.push_back(yout[2]);
        tout += 0.5;

        y[0] = yout[1];
        y[1] = yout[2];

        cout << t << ' ' << setprecision(8) << y[0] << ' ' << y[1] << endl;
    }

    cout << res[0] << ' ' << res[1] << endl;
    assert(areEqual(-11.94045, res[0]));
    assert(areEqual(3.8610102, res[1]));

    if (istate <= 0) {
        cerr << "error istate = " << istate << endl;
        throw runtime_error("Failed to compute the solution.");
    }

    return 0;
}

int test_scipy_sys(void)
{
    cout << "Running test scipy sys" << endl;
    double t, tout;
    t = 0e0;
    tout = 10;

    vector<double> y = { 10, 0 };
    int istate = 1;

    LSODA lsoda;

    // Create vector to store results. NOTE THAT yout[0] will be ignored.
    vector<double> yout;
    lsoda.lsoda_update(system_scipy, 2, y, yout, &t, tout, &istate, nullptr);

    areEqual(9.999899e+00, yout[1]);
    areEqual(-1.010111e-05, yout[2]);

    if (istate <= 0) {
        cerr << "error istate = " << istate << endl;
        exit(0);
    }
    return 0;
}

int test_fex(void)
{
    // cout << "Running test fex." << endl;
    int neq = 3;
    double t, tout;
    t = 0e0;
    tout = 0.4e0;
    vector<double> y = { 1e0, 0e0, 0.0 };
    int istate = 1;

    LSODA lsoda;
    setprecision(12);

    vector<double> res;

    vector<double> yout;
    for (size_t iout = 1; iout <= 12; iout++) {
        lsoda.lsoda_update(
            fex, neq, y, yout, &t, tout, &istate, nullptr, 1e-4, 1e-8);
        // cerr << " at t " << t << " y= " << yout[1] << ' ' << yout[2] << ' '
        // << yout[3]
        // << endl; Update the y for next iteration.
        y[0] = yout[1];
        y[1] = yout[2];
        y[2] = yout[3];

        res.push_back(y[0]);
        res.push_back(y[1]);
        res.push_back(y[2]);

        if (istate <= 0) {
            cerr << "error istate = " << istate << endl;
            exit(0);
        }
        tout = tout * 10.0E0;
    }

    vector<double> expected = { 0.985172, 3.3864e-05, 0.0147939, 0.905514,
        2.24042e-05, 0.0944634, 0.715803, 9.18446e-06, 0.284188, 0.450479,
        3.22234e-06, 0.549517, 0.183171, 8.94046e-07, 0.816828, 0.0389738,
        1.62135e-07, 0.961026, 0.00493686, 1.98442e-08, 0.995063, 0.00051665,
        2.06765e-09, 0.999483, 5.20075e-05, 2.08041e-10, 0.999948, 5.20168e-06,
        2.08068e-11, 0.999995, 5.19547e-07, 2.07819e-12, 0.999999 };

    // Assert here.
    for (size_t i = 0; i < expected.size(); i++) {
        double err = abs(expected[i] - res[i]);
        if (err > 1e-6) {
            cerr << "FAILED: Expected " << expected[i] << ". Got " << res[i]
                 << endl;
            assert(false);
        }
    }

    return 0;
}

int main(int argc, const char* argv[])
{
    (void)argc;
    (void)argv;

    auto run = [](const char* name, auto fn) {
        auto begin = chrono::steady_clock::now();
        fn();
        auto end = chrono::steady_clock::now();
        cout << "|| " << name << " time (us)= "
             << chrono::duration_cast<chrono::microseconds>(end - begin).count()
             << endl;
    };

    run("test_scipy_sys", test_scipy_sys);
    run("test_fex", test_fex);
    run("test_github_system", test_github_system);
    run("test_exponential_decay", test_exponential_decay);
    run("test_simple_harmonic_oscillator", test_simple_harmonic_oscillator);
    run("test_coupled_decay", test_coupled_decay);
    run("test_rational_ode", test_rational_ode);

    return 0;
}
