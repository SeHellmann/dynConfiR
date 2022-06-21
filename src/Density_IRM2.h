/* Density_IRM2.h - Functions for PDF calculation in the independent race model
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

#ifndef DENSITY_IRM2_H
#define DENSITY_IRM2_H

using namespace Rcpp;

#define EPSILON 1e-6

// Forward declarations

static double densIRM2 (double t, double th2, double th1,
                                    double muw, double mul, double a, double b,
                                    double wx, double wrt, double wint,
                                    double smuw, double smul,
                                    double sza, double szb, double step_width);

static double integrate_densIRM2_over_t  (double tmin, double tmax, double dt, double th2, double th1, double muw, double mul, double a, double b,
                                                     double wx, double wrt, double wint,
                                                     double smuw, double smul,
                                                     double sza, double szb, double step_width);
static double integrate_factor_b_by_szb (double b, double szb, double step_width,
                                         double t, double mul, double smul,
                                         double sigl, double tth2, double tth1);

// Main calls
NumericVector density_IRM2 (NumericVector rts, NumericVector params, int win=1, double step_width = 0.0001)
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
    double sza = params[13+win];
    double szb = params[16-win];

    muw = muw/sigmaw;
    mul = mul/sigmal;
    a = a/sigmaw;
    b = b/sigmal;
    th1 = th1/sigmal;
    th2 = th2/sigmal;
    wrt = wrt/sigmal;
    smuw = smuw/sigmaw;
    smul = smul/sigmal;
    sza = sza/sigmaw;
    szb = szb/sigmal;

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
                out[i] = fac * densIRM2( rts[i],  th2,  th1, muw,  mul,  a,  b,
                                                    wx, wrt, wint,
                                                    smuw, smul, sza, szb,
                                                    step_width);
            }
        }
    } else {
        for (int i = 0; i < length; i++) {
            if (rts[i] < 0 ) {
                out[i] = 0;
            } else {
                out[i] = fac / st0 *
                    integrate_densIRM2_over_t(rts[i]-st0, rts[i], dt,
                                                     th2,  th1, muw,  mul,  a,  b,
                                                     wx, wrt, wint, smuw, smul,
                                                     sza, szb, step_width);
            }
        }
    }
    return out;
}


static double densIRM2 (double t, double th2, double th1, double muw, double mul, double a, double b,
                                   double wx, double wrt, double wint,
                                      double smuw, double smul,
                                      double sza, double szb,
                                      double step_width) {
  double tth1 = (-th2*sqrt(t) + wrt)/(sqrt(t)*wx + wint);
  double tth2 = (-th1*sqrt(t) + wrt)/(sqrt(t)*wx + wint);
  tth2 = std::min(tth2, (double) 0);
  if (tth1 > tth2) {
    return 0;
  }

  double fac_w, fac_l1, fac_l2;
  double sigw = sqrt(2*(t+t*t*smuw*smuw));
  double sigl = sqrt(2*(t+t*t*smul*smul));
  static const double sqrtpi = 1.77245385090552;
  if (sza==0) {
    fac_w = a/(sqrtpi*sigw) * exp(-(t*muw+a)*(t*muw+a)/(sigw*sigw));
    //fac_w = a/(sqrtpi*sigw*t) * exp(-(t*muw+a)*(t*muw+a)/(sigw*sigw)); // *t
  } else {
    fac_w = muw*t*(erf((a+sza+muw*t)/sigw)-erf((a+muw*t)/sigw)) +
      sigw/sqrtpi *
      (exp(-(muw*t+a+sza)*(muw*t+a+sza)/(sigw*sigw))-
      exp(-(muw*t+a)*(muw*t+a)/(sigw*sigw)));
    // fac_w = muw/2*(erf((a+sza+muw*t)/sigw)-erf((a+muw*t)/sigw)) +   // *t
    //   sigw/(sqrtpi*2*t) *
    //   (exp(-(muw*t+a+sza)*(muw*t+a+sza)/(sigw*sigw))-
    //   exp(-(muw*t+a)*(muw*t+a)/(sigw*sigw)));
    fac_w = - fac_w / (2*sza);
  }

  //Rcpp::Rcout << "fac_w: " << fac_w <<  std::endl;
  if (szb==0) {
    fac_l1 = -erf((mul*t-tth2+b)/sigl) +
              erf((mul*t-tth1+b)/sigl);
    fac_l2 = - exp(-2*b*(mul-b*smul*smul)) *
                (erf((mul*t-tth2-b*(1+2*t*smul*smul))/sigl)-
                erf((mul*t-tth1-b*(1+2*t*smul*smul))/sigl));
    // fac_l = 1 - erf((b+mul*t)/sigl) -
    //   exp(2*b*(smul*smul*b-mul)) * (1+erf((b-t*mul+2*t*b*smul*smul)/sigl));
  } else {
    fac_l1 = - (   (b+szb+mul*t-tth2)*erf((b+szb+mul*t-tth2)/sigl) - (b+mul*t-tth2)*erf((b+mul*t-tth2)/sigl) +
      sigl/(sqrtpi) * (exp(-(b+szb+mul*t-tth2)*(b+szb+mul*t-tth2)/(sigl*sigl)) -
      exp(-(b+mul*t-tth2)*(b+mul*t-tth2)/(sigl*sigl))) ) +
      (   (b+szb+mul*t-tth1)*erf((b+szb+mul*t-tth1)/sigl) - (b+mul*t-tth1)*erf((b+mul*t-tth1)/sigl) +
      sigl/(sqrtpi) * (exp(-(b+szb+mul*t-tth1)*(b+szb+mul*t-tth1)/(sigl*sigl)) -
      exp(-(b+mul*t-tth1)*(b+mul*t-tth1)/(sigl*sigl)))   ) ;

    // fac_l = szb - ((b+szb+mul*t)*erf((b+szb+mul*t)/sigl) - (b+mul*t)*erf((b+mul*t)/sigl) +
    //   sigl/(sqrtpi) * (exp(-(b+szb+mul*t)*(b+szb+mul*t)/(sigl*sigl)) -
    //   exp(-(b+mul*t)*(b+mul*t)/(sigl*sigl)))  );
    if (smul==0) {
      fac_l2 = 1/(2*mul)*(
         ( exp(2*mul*tth2)*(erf((b+szb+mul*t+tth2)/sigl)-erf((b+mul*t+tth2)/sigl)) -
        exp(-2*(b+szb)*mul)*erf((b+szb-mul*t+tth2)/sigl) +
        exp(-2*(b)*mul)*erf((b-mul*t+tth2)/sigl) ) -
        (   exp(2*mul*tth1)*(erf((b+szb+mul*t+tth1)/sigl)-erf((b+mul*t+tth1)/sigl)) -
        exp(-2*(b+szb)*mul)*erf((b+szb-mul*t+tth1)/sigl) +
        exp(-2*(b)*mul)*erf((b-mul*t+tth1)/sigl) )  ) ;
      // fac_l = fac_l - 1/(2*mul)* ( erf((b+szb+mul*t)/sigl)-erf((b+mul*t)/sigl) -
      //   exp(-2*(b+szb)*mul)*(erf((b+szb-mul*t)/sigl)+1) +
      //   exp(-2*(b)*mul)*(erf((b-mul*t)/sigl)+1) );
    } else {
      fac_l2 = integrate_factor_b_by_szb(b,szb, step_width,
                                                t, mul, smul, sigl, tth2, tth1);
    }
    fac_l1 = fac_l1 / szb;
    fac_l2 = fac_l2 / szb;

  }
  //Rcpp::Rcout << "fac_l: " << fac_l <<  std::endl;
  return -fac_w*(fac_l1-fac_l2)/(2*t);
}

static double integrate_densIRM2_over_t (double tmin, double tmax, double dt, double th2, double th1, double muw, double mul, double a, double b,
                                                    double wx, double wrt, double wint,
                                                       double smuw, double smul,
                                                       double sza, double szb,
                                                       double step_width) {
    double x;
    double result = 0;
    for(x = tmin+0.5*dt; x < tmax; x += dt)
    {
        if ( x > 0)
        {
        result +=  dt * densIRM2 ( x,  th2,  th1, muw,  mul,  a,  b, wx, wrt, wint, smuw, smul, sza, szb,
                                   step_width);
        }
    }
    return result;
}

static double integrate_factor_b_by_szb (double b, double szb, double step_width,
                                         double t, double mul, double smul,
                                         double sigl, double tth2, double tth1) {
  int nsteps;
  double db;
  nsteps = std::max(4, (int) (szb/step_width));
  db = szb / nsteps;

  double B;
  double result = 0;
  for(B = b+0.5*db; B < b+szb; B += db)
  {
    result +=  db*(exp(-2*B*(mul-B*smul*smul))*
      (erf(((2*t*smul*smul+1)*B -t*mul+tth2)/sigl) -
      erf(((2*t*smul*smul+1)*B -t*mul+tth1)/sigl)));
  }
  return result;
}

#endif // DENSITY_IRM2_H
