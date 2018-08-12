/***
 *       Filename:  LSODE.h
 *
 *    Description:  See LSODE.cpp file to more information.
 *
 *        Version:  0.0.1
 *        Created:  2018-08-09

 *       Revision:  none
 *
 *         Author:  Dilawar Singh <dilawars@ncbs.res.in>
 *   Organization:  NCBS Bangalore
 *
 *        License:  GNU GPL3
 */

#ifndef LSODE_H
#define LSODE_H

#include <memory>
#include <array>
#include <cmath>
#include <vector>
#include <memory>

using namespace std;

/* --------------------------------------------------------------------------*/
/**
 * @Synopsis  Type definition of LSODA ode system. See the file test_LSODA.cpp
 * for an example.
 *
 * @Param time, double
 * @Param y, array of double.
 * @Param dydt, array of double
 * @Param data, void*
 *
 * @Returns void
 */
/* ----------------------------------------------------------------------------*/
typedef void  (*LSODA_ODE_SYSTEM_TYPE) (double t, double * y, double * dydt, void *);


class LSODA
{

public:

    LSODA( );
    ~LSODA();

    int idamax( const size_t n, const double* const, const size_t incx);

    void dscal(const size_t n, const double da, double* dx, const size_t incx);

    double ddot(const size_t n, const double* const dx, const size_t ncx
            , const double* const dy, const size_t incy
            );

    void daxpy(const size_t n, const double da, const double* const dx
            , const size_t incx, double* dy, const size_t incy
            );

    void dgesl(double** const a, const size_t n, const int* const ipvt, double* b
            , const size_t job
            );

    void dgefa(double** const a, const size_t n, int* const ipvt, size_t* const info);

    void prja(const size_t neq, double* const y, LSODA_ODE_SYSTEM_TYPE f, void *_data);

    int n_lsoda(double y[], int n, double *x, double xout, double eps
            , const double yscal[], LSODA_ODE_SYSTEM_TYPE devis, void *data
            );

    void lsoda( LSODA_ODE_SYSTEM_TYPE f, int neq
                , double *y, double *t, double tout
                , int itask, int *istate, int iopt, int jt
                , int iwork1, int iwork2, int iwork5, int iwork6, int iwork7, int iwork8, int iwork9
                , double rwork1, double rwork5, double rwork6, double rwork7
                , void *_data
              );

    void correction( const size_t neq, double* const y
            , LSODA_ODE_SYSTEM_TYPE f, int *corflag
            , double pnorm, double *del, double *delp, double *told
            , int *ncf, double *rh, int *m, void *_data
            );

    void stoda(int neq, double *y, LSODA_ODE_SYSTEM_TYPE f, void *_data);

    // We call this function in VoxelPools::
    void lsoda_update( LSODA_ODE_SYSTEM_TYPE f, const size_t neq
            , const double* const y, std::vector<double>& yout
            , double* t, const double tout
            , int* istate
            , void* const _data
            , double rtol = 1e-4, double atol = 1e-6   // Tolerance
            );

    void     terminate(int *istate);
    void     terminate2(double *y, double *t);
    void     successreturn(double *y, double *t, int itask, int ihit, double tcrit, int *istate);
    void     _freevectors(void);
    void     ewset(double *ycur);
    void     resetcoeff(void);
    void     solsy(double *y);
    void     endstoda(void);
    void     orderswitch(double *rhup, double dsm, double *pdh, double *rh, int *orderflag);
    void     intdy(double t, int k, double *dky, int *iflag);
    void     corfailure(double *told, double *rh, int *ncf, int *corflag);
    void     methodswitch(double dsm, double pnorm, double *pdh, double *rh);
    void     cfode(int meth);
    void     scaleh(double *rh, double *pdh);
    double   fnorm(int n, double **a, double *w);
    double   vmnorm(int n, double *v, double *w);

private:
    int      g_nyh = 0, g_lenyh = 0;

    int      ml, mu, imxer;
    double   sqrteta, *yp1, *yp2;

    // NOTE: initialize in default constructor. Older compiler e.g. 4.8.4 would
    // produce error if these are initialized here.

    array<int,3> mord; // = {0, 12, 5};

    // = {  0., 0.5,  0.575, 0.55, 0.45, 0.35, 0.25,  0.2, 0.15, 0.1,  0.075, 0.05, 0.025 };
    array<double, 13>  sm1;

    array<double, 14> el; // = {0};
    array<double, 13> cm1; // = {0};
    array<double,6> cm2; // = {0};

    array<array<double, 14>, 13> elco;
    array<array<double,4>, 13> tesco;

    int      illin = 0, init = 0, mxstep, mxhnil, nhnil, ntrep = 0, nslast, nyh, ierpj, iersl,
             jcur, jstart, kflag, l, meth, miter, maxord, maxcor, msbp, mxncf, n, nq, nst,
             nfe, nje, nqu;
    int      ixpr = 0, jtyp, mused, mxordn, mxords;

    double   ccmax, el0, h, hmin, hmxi, hu, rc, tn;
    double   tsw, pdnorm;
    double   conit, crate, hold, rmax;

    int      ialth, ipup, lmax, nslp;
    double   pdest, pdlast, ratio;
    int      icount, irflag;

    double **yh, **wm, *ewt, *savf, *acor;
    int     *ipvt;

private:
    int itol_ = 2;
    int istate_ = 1;
    std::vector<double> rtol_;
    std::vector<double> atol_;

public:
    void* param = nullptr;
};


#endif /* end of include guard: LSODE_H */
