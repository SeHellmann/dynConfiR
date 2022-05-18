/* Density_IRM.h - Functions for PDF calculation in the independent race model
 *
 * Copyright (C) 2022 Sebastian Hellmann.
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation; either version 3 of the
 * License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 */

#ifndef DENSITY_IRM_H
#define DENSITY_IRM_H

using namespace Rcpp;

#define EPSILON 1e-6

// Forward declarations

static double densIRM_equalbounds (double t, double th2, double th1, double muw, double mul, double a, double b, double wx, double wrt, double wint);
static double densIRM_differbounds (double t, double th2, double th1,
                             double muw, double mul,
                             double wx, double wrt, double wint,
                             NumericVector C, NumericVector expC, NumericVector Xis);

static double integrate_densIRM_equalbounds_over_t  (double tmin, double tmax, double dt, double th2, double th1, double muw, double mul, double a, double b,
                                                     double wx, double wrt, double wint);
static double integrate_densIRM_differbounds_over_t (double tmin, double tmax, double dt, double th2, double th1, double muw, double mul,
                                                     double wx, double wrt, double wint,
                                                     NumericVector C, NumericVector expC, NumericVector Xis);


// Main calls
NumericVector density_IRM (NumericVector rts, NumericVector params, int win=1, double step_width = 0.0001)
{
    int length = rts.length();
    NumericVector out(length, 0.0);

    double muw = params[win-1];
    double mul = params[2-win];
    double a = params[1+win];
    double b = params[4-win];
    double sigmaw = params[3+win];
    double sigmal = params[6-win];
    double st0 = params[8];
    double th1 = params[6];
    double th2 = params[7];
    double wx = params[9];
    double wrt = params[10];
    double wint = params[11];

    muw = muw/sigmaw;
    mul = mul/sigmal;
    a = a/sigmaw;
    b = b/sigmal;
    th1 = th1/sigmal;
    th2 = th2/sigmal;
    wrt = wrt/sigmal;

    int nsteps;
    double dt;
    if (st0 < EPSILON) {
        st0 = 0;
    }
    if (st0 == 0) {
        nsteps=1;
        dt=0;
    } else {
        nsteps = std::max(4, (int) (st0/step_width));
        dt = st0 / nsteps;
    }

    static const float inv_sqrt_2pi = 0.3989422804014327;

    if (a==b) {
        double fac = - inv_sqrt_2pi*a/2;
        if (st0==0) {
            for (int i = 0; i < length; i++) {
                if (rts[i] < 0 ) {
                    out[i] = 0;
                } else {
                    out[i] = fac * densIRM_equalbounds( rts[i],  th2,  th1, muw,  mul,  a,  b,
                                                        wx, wrt, wint);
                }
            }
        } else {
            for (int i = 0; i < length; i++) {
                if (rts[i] < 0 ) {
                    out[i] = 0;
                } else {
                    out[i] = fac / st0 *
                        integrate_densIRM_equalbounds_over_t(rts[i]-st0, rts[i], dt,
                                                         th2,  th1, muw,  mul,  a,  b,
                                                         wx, wrt, wint);
                }
            }
        }
    } else {
        double fac = inv_sqrt_2pi/4;
        NumericVector C = NumericVector::create(1, -1, -1, 1);
        NumericVector expC1 = NumericVector::create(a, 0, a);
        NumericVector expC2 = NumericVector::create(0, b, b);

        NumericVector expC = -2 * (muw*expC1 + mul*expC2);

        expC.push_front(0);

        NumericVector Xis = NumericVector::create(a, -a, a, -a,
                             b, b, -b, -b);
        Xis.attr("dim") = Dimension(4,2);
        if (st0==0) {
            for (int i = 0; i < length; i++) {
                if (rts[i] < 0 ) {
                    out[i] = 0;
                } else {
                    out[i] = fac *  densIRM_differbounds (rts[i], th2,  th1, muw,  mul, wx, wrt, wint,
                                       C, expC,  Xis);
                }
            }
        } else {
            for (int i = 0; i < length; i++) {
                if (rts[i] < 0 ) {
                    out[i] = 0;
                } else {
                    out[i] = fac / st0 *
                    integrate_densIRM_differbounds_over_t(rts[i]-st0, rts[i], dt,
                                                         th2,  th1, muw,  mul, wx, wrt, wint,
                                                         C, expC, Xis);
                }
            }
        }
    }
    return out;
}


static double densIRM_equalbounds (double t, double th2, double th1, double muw, double mul, double a, double b,
                                   double wx, double wrt, double wint) {
    double sig2t = 2*t;
    double tth1 = (-th2*sqrt(t) + wrt)/(sqrt(t)*wx + wint);
    double tth2 = (-th1*sqrt(t) + wrt)/(sqrt(t)*wx + wint);
    tth2 = std::min(tth2, (double) 0);
    if (tth1 > tth2) {
        return 0;
    }
    double res = 1/pow(t, 1.5) *
        exp(-(muw*t+a)*(muw*t+a)/sig2t) *
        (erf((tth2-(mul*t+a))/sqrt(sig2t))-erf((tth1-(mul*t+a))/sqrt(sig2t)) -
        exp(-2*b*mul)*(erf((tth2-(mul*t-a))/sqrt(sig2t))-erf((tth1-(mul*t-a))/sqrt(sig2t))));
    return res;
}

static double integrate_densIRM_equalbounds_over_t (double tmin, double tmax, double dt, double th2, double th1, double muw, double mul, double a, double b,
                                                    double wx, double wrt, double wint) {
    double x;
    double result = 0;
    for(x = tmin+0.5*dt; x < tmax; x += dt)
    {
        if ( x > 0)
        {
        result +=  dt * densIRM_equalbounds ( x,  th2,  th1, muw,  mul,  a,  b, wx, wrt, wint);
        }
    }
    return result;
}

static double densIRM_differbounds (double t, double th2, double th1,
                             double muw, double mul, double wx, double wrt, double wint,
                             NumericVector C, NumericVector expC, NumericVector Xis){
    double sig2t = 2*t;
    double tth1 = (-th2*sqrt(t) + wrt)/(sqrt(t)*wx + wint);
    double tth2 = (-th1*sqrt(t) + wrt)/(sqrt(t)*wx + wint);
    tth2 = std::min(tth2, (double) 0);
    if (tth1 > tth2) {
        return 0;
    }
    double temp = 0;
    for (int j = 0; j < 4; j++) {
        double x1tilde = -Xis(j,0)-muw*t;
        temp = temp + C[j]* exp(expC[j] - (x1tilde)*(x1tilde)/(sig2t)) * x1tilde *
            (erf((tth2-Xis(j,1)-mul*t)/sqrt(sig2t)) - erf((tth1-Xis(j,1)-mul*t)/sqrt(sig2t)));
    }
    double res = temp/pow(t, 1.5);
    return res;
}

static double integrate_densIRM_differbounds_over_t (double tmin, double tmax, double dt, double th2, double th1, double muw, double mul,
                                                     double wx, double wrt, double wint,
                                                     NumericVector C, NumericVector expC, NumericVector Xis) {
    double x;
    double result = 0;
    for(x = tmin+0.5*dt; x < tmax; x += dt)
    {
        if ( x > 0)
        {
        result +=  dt * densIRM_differbounds ( x,  th2,  th1, muw,  mul,  wx, wrt, wint, C, expC, Xis);
        }
    }
    return result;
}

#endif // DENSITY_IRM_H
