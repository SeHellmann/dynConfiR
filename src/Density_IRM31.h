/* Density_IRM3.h - Functions for PDF calculation in the independent race model
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

#ifndef DENSITY_IRM3_H
#define DENSITY_IRM3_H

using namespace Rcpp;

#define EPSILON 1e-6

// Forward declarations

static double densIRM3 (double t, double th2, double th1,
                        double a, double b,
                        double muw, double mul,
                        double wx, double wrt, double wint,
                        double smuw, double smul);

static double integrate_densIRM3_over_t (double tmin, double tmax,
                                         double dt, double th2, double th1,
                                         double a, double b,
                                         double muw, double mul,
                                         double wx, double wrt, double wint,
                                         double smuw, double smul);

// Main calls
NumericVector density_IRM3 (NumericVector rts, NumericVector params, int win=1, double step_width = 0.0001)
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
    double smuw = params[11+win];
    double smul = params[14-win];


    muw = muw/sigmaw;
    mul = mul/sigmal;
    a = a/sigmaw;
    b = b/sigmal;
    th1 = th1/sigmal;
    th2 = th2/sigmal;
    wrt = wrt/sigmal;
    smuw = smuw/sigmaw;
    smul = smul/sigmal;


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
    double fac = - inv_sqrt_2pi*a/2;

    if (st0==0) {
      for (int i = 0; i < length; i++) {
        if (rts[i] < 0 ) {
          out[i] = 0;
        } else {
          out[i] = fac *   // 1/sqrt(2*pi)
            densIRM3 (rts[i], th2,  th1, a, b, muw,  mul, wx, wrt, wint,
                      smuw, smul);
        }
      }
    } else {
      for (int i = 0; i < length; i++) {
        if (rts[i] < 0 ) {
          out[i] = 0;
        } else {
          out[i] = fac / st0 *    // 1/sqrt(2*pi)
            integrate_densIRM3_over_t(rts[i]-st0, rts[i], dt,
                                      th2,  th1, a, b, muw,  mul, wx, wrt, wint,
                                      smuw, smul);
        }
      }
    }

    return out;
}


static double densIRM3 (double t, double th2, double th1, double a, double b,
                             double muw, double mul, double wx, double wrt, double wint,
                             double smuw, double smul
                             ){
  //double facA, facB, res;
  double Sig2tw = 2.0*(t + t*t*smuw*smuw);
  double Sig2tl = 2.0*(t + t*t*smul*smul); //sqrt(2.0*(t + t*t*smul *smul));

  double tth1 = (-th2*sqrt(t) + wrt)/(sqrt(t)*wx + wint);
  double tth2 = (-th1*sqrt(t) + wrt)/(sqrt(t)*wx + wint);
  tth2 = std::min(tth2, (double) 0);
  if (tth1 > tth2) {
    return 0; // ensure tth1 < tth2 < 0 (i.e. tth represent distance x2 from 0)
  }
  // facA = -a/(2.0*t*sqrt(Sigw2))*
  //     exp(-(a+t*muw)*(a+t*muw)/(2*Sigw2));
  // // facB = 2.0*((
  // //   Phi((tth2-t*mul+b)/sqrt(Sigl2*t))-
  // //     Phi((tth1-t*mul+b)/sqrt(Sigl2*t))) +
  // //     exp(2*b*(smul*smul*b - mul))*
  // //     (Phi((tth2-t*mul+b*(1+2*smul*smul*t))/sqrt(Sigl2*t)) -
  // //     Phi((tth1-t*mul+b*(1+2*smul*smul*t))/sqrt(Sigl2*t))));
  // facB = (
  //   erf((tth2-t*mul-b)/Sigl)-
  //     erf((tth1-t*mul-b)/Sigl)) +
  //     exp(2*b*(smul*smul*b - mul))*
  //     (erf((tth2-t*mul+b*(1+2*smul*smul*t))/Sigl) -
  //     erf((tth1-t*mul+b*(1+2*smul*smul*t))/Sigl));
  //
  // res = facA*facB;

  // a/(2*t) * 1/(sqrt(2*pi*(t+t^2*smu1^2))) * exp(-(-a+t*mu1)^2/(2*(t+t^2*smu1^2))) *
  //   (erf((-th1 - (t*mu2-b))/(sqrt(2*(t+t*t*smu2*smu2)))) -
  //   erf((-th2 - (t*mu2-b))/(sqrt(2*(t+t*t*smu2*smu2)))) -
  //   exp(2*b*(mu2+b*smu2^2)) *
  //   (erf((-th1 - t*mu2-b*(1+2*t*smu2^2))/(sqrt(2*(t+t*t*smu2*smu2))))-
  //   erf((-th2 - t*mu2-b*(1+2*t*smu2^2))/(sqrt(2*(t+t*t*smu2*smu2))))))
  //
  double res = 1/(t*sqrt(t+t*t*smuw*smuw)) *
    exp(-(muw*t+a)*(muw*t+a)/Sig2tw) *
    (erf((tth2-(mul*t+b))/sqrt(Sig2tl))-erf((tth1-(mul*t+b))/sqrt(Sig2tl)) -
    exp(-2*b*mul + 2*b*b*smul*smul +
      log(erf((tth2-(mul*t-b*(1+2*smul*smul*t)))/sqrt(Sig2tl)) -
      erf((tth1-(mul*t-b*(1+2*smul*smul*t)))/sqrt(Sig2tl)))));
  return res;
}

static double integrate_densIRM3_over_t (double tmin, double tmax, double dt, double th2, double th1,
                                         double a, double b,
                                         double muw, double mul,
                                         double wx, double wrt, double wint,
                                         double smuw, double smul) {
  double x;
  double result = 0;
  for(x = tmin+0.5*dt; x < tmax; x += dt)
  {
    if ( x > 0)
    {
      result +=  dt * densIRM3 ( x,  th2,  th1, a, b, muw,  mul,  wx, wrt, wint, smuw, smul);
    }
  }
  return result;
}

#endif // DENSITY_IRM3_H
