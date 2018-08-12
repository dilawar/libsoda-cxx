/***
 *       Filename:  LSODE.h_
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
typedef void  (*LSODA_ODE_SYSTEM_TYPE) (double t, double* y, double* dydt, void *);


class LSODA
{

public:

    LSODA( );
    ~LSODA();

    int idamax( const size_t n, const double* const, const size_t incx);

    void dscal(const size_t n, const double da, double* dx, const size_t incx);

    double ddot(const size_t n, const double* const dx, const int ncx
            , const double* const dy, const int incy
            );

    void daxpy(const size_t n, const double da, const double* const dx
            , const int incx, double* dy, const int incy
            );

    void dgesl( const vector<vector<double>>& a, const size_t n, vector<int>& ipvt, double* b
            , const size_t job
            );

    void dgefa( vector<vector<double>>& a, const size_t n, vector<int>& ipvt, size_t* const info);

    void prja(const size_t neq, vector<double>& y, LSODA_ODE_SYSTEM_TYPE f, void *_data);

    void lsoda( LSODA_ODE_SYSTEM_TYPE f, const size_t neq
                , vector<double>& y, double *t, double tout
                , int itask, int *istate, int iopt, int jt
                , array<int, 7>& iworks, array<double, 4>& rworks
                , void *_data
              );

    void correction( const size_t neq, vector<double>& y
            , LSODA_ODE_SYSTEM_TYPE f, size_t *corflag
            , double pnorm
            , double *del, double *delp, double *told
            , size_t *ncf, double *rh, size_t *m, void *_data
            );

    void stoda(const size_t neq, vector<double>& y, LSODA_ODE_SYSTEM_TYPE f, void *_data);

    // We call this function in VoxelPools::
    void lsoda_update( LSODA_ODE_SYSTEM_TYPE f, const size_t neq
            , vector<double>& y, std::vector<double>& yout
            , double* t, const double tout
            , int* istate
            , void* const _data
            , double rtol = 1e-4, double atol = 1e-6   // Tolerance
            );

    void     terminate(int *istate);
    void     terminate2( vector<double>& y, double *t);
    void     successreturn( vector<double>& y, double *t, int itask, int ihit, double tcrit, int *istate);
    void     _freevectors(void);
    void     ewset( const vector<double>& ycur);
    void     resetcoeff(void);
    void     solsy(double *y);
    void     endstoda(void);
    void     orderswitch(double *rhup, double dsm, double *pdh, double *rh, size_t *orderflag);
    void     intdy(double t, int k, vector<double>& dky, int *iflag);
    void     corfailure(double *told, double *rh, size_t *ncf, size_t *corflag);
    void     methodswitch(double dsm, double pnorm, double *pdh, double *rh);
    void     cfode(int meth_);
    void     scaleh(double *rh, double *pdh);
    double   fnorm(int n, const vector<vector<double>>& a, const vector<double>& w);
    double   vmnorm(const size_t n, const vector<double>& v, const vector<double>& w);

private:

    size_t   ml, mu, imxer;
    double   sqrteta;

    // NOTE: initialize in default constructor. Older compiler e.g. 4.8.4 would
    // produce error if these are initialized here. With newer compiler,
    // initialization can be done here.
    array<size_t,3> mord;
    array<double, 13>  sm1;

    array<double, 14> el; // = {0};
    array<double, 13> cm1; // = {0};
    array<double,6> cm2; // = {0};

    array<array<double, 14>, 13> elco;
    array<array<double,4>, 13> tesco;

    size_t  illin, init, ierpj, iersl, jcur, l, miter, maxord, maxcor, msbp, mxncf;

    int      kflag, jstart;

    size_t   ixpr = 0, jtyp, mused, mxordn, mxords = 12;
    size_t  meth_;

    size_t   n, nq, nst, nfe, nje, nqu;
    size_t   mxstep, mxhnil;
    size_t   nslast, nhnil, ntrep, nyh;

    double   ccmax, el0, h_=.0;
    double   hmin, hmxi, hu, rc, tn_=0.0;
    double   tsw, pdnorm;
    double   conit, crate, hold, rmax;

    size_t   ialth, ipup, lmax;
    size_t   nslp;
    double   pdest, pdlast, ratio;
    int      icount, irflag;

    vector<double> ewt;
    vector<double> savf;
    vector<double> acor;
    vector<vector<double>> yh_;
    vector<vector<double>> wm_;

    vector<int> ipvt;

private:
    int itol_ = 2;
    std::vector<double> rtol_;
    std::vector<double> atol_;

public:
    void* param = nullptr;
};


#endif /* end of include guard: LSODE_H */
