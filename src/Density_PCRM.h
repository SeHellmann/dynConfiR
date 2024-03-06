/* Density_PCRM.h - Functions for PDF calculation in the partially anti-correlated race model
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


#ifndef DENSITY_PCRM_H
#define DENSITY_PCRM_H

using namespace Rcpp;

#define EPSILON 1e-6

// Forward declarations

static double densPCRM(double t, double th2, double th1,
                             double muw, double mul, double wx, double wrt, double wint,
                             NumericVector C, NumericVector expC, NumericVector Xis);
static double integrate_densPCRM_over_t (double tmin, double tmax, double dt, double th2, double th1, double muw, double mul,
                                         double wx, double wrt, double wint,
                                         NumericVector C, NumericVector expC, NumericVector Xis);

double fac_errf = 2.170803763674803; // sqrt(M_PI*3/2); used in densPCRM

// Main calls
NumericVector density_PCRM (NumericVector rts, NumericVector params, int win=1, double step_width = 0.0001)
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

    double fac = 1/(sqrt(3.)*4*M_PI);
    NumericVector C = NumericVector::create(1, -1, -1, 1, 1, -1);
    NumericVector expC1 = NumericVector::create(0, a, 0, a, a+b, a+b);
    NumericVector expC2 = NumericVector::create(0, 0, b, a+b, b, a+b);

    NumericVector expC = -2* (muw*expC1 + mul*expC2);

    NumericVector Xis = NumericVector::create(a, -a, a+b, b, -a-b, -b,
                         b, a+b, -b, -a-b, a, -a);
    Xis.attr("dim") = Dimension(6,2);

    if (st0==0) {

        for (int i = 0; i < length; i++) {
            if (rts[i] < 0 ) {
                out[i] = 0;
            } else {
                out[i] = fac *  densPCRM (rts[i], th2,  th1, muw,  mul, wx, wrt, wint,
                                   C, expC,  Xis);
              if (i % 200 ==0 ) Rcpp::checkUserInterrupt();
            }
        }
    } else {
        for (int i = 0; i < length; i++) {
            if (rts[i] < 0 ) {
                out[i] = 0;
            } else {
                out[i] = fac / st0 *
                    integrate_densPCRM_over_t(rts[i]-st0, rts[i], dt,
                                              th2,  th1, muw,  mul, wx, wrt, wint,
                                              C, expC, Xis);
              if (i % 200 ==0 ) Rcpp::checkUserInterrupt();
            }
        }
    }

    return out;
}

static double densPCRM (double t, double th2, double th1,
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
    for (int j = 0; j < 6; j++) {
        double num_erfn2 = tth2 - Xis(j,1) - mul*t - (Xis(j,0)+muw*t)/2;
        double num_erfn1 = tth1 - Xis(j,1) - mul*t - (Xis(j,0)+muw*t)/2;
        double den_phi = sig2t * 0.75;
        temp = temp + C[j]* exp(expC[j] - (Xis(j, 0)+muw*t)*(Xis(j, 0)+muw*t)/(sig2t)) *
            (fac_errf/sqrt(t)* (- (Xis(j,0)+muw*t)) * (erf(num_erfn2/sqrt(den_phi)) - erf(num_erfn1/sqrt(den_phi))) -
            (exp(-num_erfn2*num_erfn2/den_phi) - exp(-num_erfn1*num_erfn1/den_phi)));
    }

    double res = temp/t;
    return res;
}

static double integrate_densPCRM_over_t (double tmin, double tmax, double dt, double th2, double th1, double muw, double mul,
                                         double wx, double wrt, double wint,
                                         NumericVector C, NumericVector expC, NumericVector Xis) {
    double x;
    double result = 0;
    for(x = tmin+0.5*dt; x < tmax; x += dt)
    {
        if ( x > 0)
        {
        result += dt * densPCRM ( x,  th2,  th1, muw, mul, wx, wrt, wint, C, expC, Xis);
        }
    }
    return  result;
}

#endif // DENSITY_PCRM_H
