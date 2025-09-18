/* SeqSampConf.cpp - Main source file for the RCpp implementation of the
 * sequential sampling models contained in the dynConfiR package
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


#include <Rcpp.h>
//#include <R_ext/Rdynload.h>
#include <iostream>
#include <sstream>

#include "SeqSampConf.h"

using namespace Rcpp;

// R-callable PDF for 2DSD - pass boundary to retrieve (1 = lower, 2 = upper)
// [[Rcpp::export]]
NumericVector d_2DSD (NumericVector rts, NumericVector params, double precision=1e-5, int boundary=2, bool stop_on_error=true, int stop_on_zero=false)
{
    int length = rts.length();
    if (length > MAX_INPUT_VALUES) { Rcpp::stop("Number of RT values passed in exceeds maximum of %d.\n", MAX_INPUT_VALUES); }

    if ((boundary < 1) || (boundary > 2)) { Rcpp::stop ("Boundary must be either 2 (upper) or 1 (lower)\n"); }

    NumericVector out(length, 0.0);  // Should default to 0s when creating NumericVector, but just in case..

    if (!ValidateParams(params, true))
    {
        if (stop_on_error) { Rcpp::stop("Error validating parameters.\n"); }
                      else { return out; }
    }

    // Add tuning values for numerical integrations at the end of parameters
    // ToDo: Optimize and check precision values
    if (precision >= 1) {
      params.push_back(0.0089045 * exp(-1.037580*(precision-3.5))); // TUNE_INT_T0
      params.push_back(0.0508061 * exp(-1.022373*(precision-3.5))); // TUNE_INT_Z
      //     These have been added to optimise code paths by treating very small variances as 0
      //     e.g. with precision = 3, sv or sz values < 10^-5 are considered 0
      params.push_back(pow (10, -(precision+2.0))); // TUNE_SZ_EPSILON
      params.push_back(pow (10, -(precision+2.0))); // TUNE_ST0_EPSILON
    } else {
      params.push_back(precision); // TUNE_INT_T0
      params.push_back(precision); // TUNE_INT_Z
      params.push_back(0); // TUNE_SZ_EPSILON
      params.push_back(0); // TUNE_ST0_EPSILON
    }

    out = density_2DSD (rts, params, boundary-1, stop_on_zero);

    //delete g_Params;
    return out;
}


// [[Rcpp::export]]
NumericVector d_WEVmu (NumericVector rts, NumericVector params, double precision=1e-5, int boundary=2,
                       bool stop_on_error=true, int stop_on_zero=false)
{
    int length = rts.length();
    if (length > MAX_INPUT_VALUES) { Rcpp::stop("Number of RT values passed in exceeds maximum of %d.\n", MAX_INPUT_VALUES); }

    if ((boundary < 1) || (boundary > 2)) { Rcpp::stop ("Boundary must be either 2 (upper) or 1 (lower)\n"); }

    NumericVector out(length, 0.0);  // Should default to 0s when creating NumericVector, but just in case..

    if (!ValidateParams(params, true))
    {
      if (stop_on_error) { Rcpp::stop("Error validating parameters.\n"); }
      else { return out; }
    }

    if (precision >= 1) {
      params.push_back(0.0089045 * exp(-1.037580*(precision-3.5))); // TUNE_INT_T0
      params.push_back(0.0508061 * exp(-1.022373*(precision-3.5))); // TUNE_INT_Z
      //     These have been added to optimise code paths by treating very small variances as 0
      //     e.g. with precision = 3, sv or sz values < 10^-5 are considered 0
      params.push_back(pow (10, -(precision+2.0))); // TUNE_SZ_EPSILON
      params.push_back(pow (10, -(precision+2.0))); // TUNE_ST0_EPSILON
    } else {
      params.push_back(precision); // TUNE_INT_T0
      params.push_back(precision); // TUNE_INT_Z
      params.push_back(0); // TUNE_SZ_EPSILON
      params.push_back(0); // TUNE_ST0_EPSILON
    }

    out = density_WEVmu (rts, params, boundary-1, stop_on_zero);

    return out;
}

// R-callable PDF for DDConf - pass boundary to retrieve (1 = lower, 2 = upper)
// [[Rcpp::export]]
NumericVector d_DDConf (NumericVector rts, NumericVector params, double precision=6, int boundary=2,
                         bool stop_on_error=true, bool stop_on_zero=false,
                         double st0precision=0.01)
{
  int length = rts.length();
  if (length > MAX_INPUT_VALUES) { Rcpp::stop("Number of RT values passed in exceeds maximum of %d.\n", MAX_INPUT_VALUES); }

  if ((boundary < 1) || (boundary > 2)) { Rcpp::stop ("Boundary must be either 2 (upper) or 1 (lower)\n"); }


  NumericVector out(length, 0.0);  // Should default to 0s when creating NumericVector, but just in case..

  if (!ValidateParams(params, true))
  {
    if (stop_on_error) { Rcpp::stop("Error validating parameters.\n"); }
    else { return out; }
  }

  // Add tuning values for numerical integrations at the end of parameters
  // ToDo: Optimize and check precision values
  params.push_back(0.0089045 * exp(-1.037580*precision)); // TUNE_INT_T0
  params.push_back(0.0508061 * exp(-1.022373*precision)); // TUNE_INT_Z
  //     These have been added to optimise code paths by treating very small variances as 0
  //     e.g. with precision = 3, sv or sz values < 10^-5 are considered 0
  params.push_back(pow (10, -(precision+2.0))); // TUNE_SZ_EPSILON
  params.push_back(pow (10, -(precision+2.0))); // TUNE_ST0_EPSILON

  out = density_DDConf (rts, params, boundary-1, stop_on_zero, st0precision);

  return out;
}

// [[Rcpp::export]]
NumericVector d_IRM (NumericVector rts, NumericVector params, int win=1,  double step_width = 0.0001)
{
  int length = rts.length();

  if (params.length()!= 12) {Rcpp::stop ("Wrong number of parameters given. (Must be 12)\n"); }
  if ((win < 1) || (win > 2)) { Rcpp::stop ("Boundary must be either 2 (upper) or 1 (lower)\n"); }

  NumericVector out(length, 0.0);  // Should default to 0s when creating NumericVector, but just in case..

  out = density_IRM (rts, params, win, step_width);

  return out;
}

// [[Rcpp::export]]
NumericVector d_IRM2 (NumericVector rts, NumericVector params, int win=1,  double step_width = 0.0001)
{
  int length = rts.length();

  if (params.length()!= 16) {Rcpp::stop ("Wrong number of parameters given. (Must be 16)\n"); }
  if ((win < 1) || (win > 2)) { Rcpp::stop ("Boundary must be either 2 (upper) or 1 (lower)\n"); }

  NumericVector out(length, 0.0);  // Should default to 0s when creating NumericVector, but just in case..

  out = density_IRM2 (rts, params, win, step_width);

  return out;
}


// [[Rcpp::export]]
NumericVector d_IRM3 (NumericVector rts, NumericVector params, int win=1,  double step_width = 0.0001)
{
  int length = rts.length();

  if (params.length()!= 14) {Rcpp::stop ("Wrong number of parameters given. (Must be 14)\n"); }
  if ((win < 1) || (win > 2)) { Rcpp::stop ("Boundary must be either 2 (upper) or 1 (lower)\n"); }

  NumericVector out(length, 0.0);  // Should default to 0s when creating NumericVector, but just in case..

  out = density_IRM3 (rts, params, win, step_width);

  return out;
}



// [[Rcpp::export]]
NumericVector d_PCRM (NumericVector rts, NumericVector params, int win=1, double step_width = 0.0001)
{
    int length = rts.length();

    if (params.length()!= 12) {Rcpp::stop ("Wrong number of parameters given. (Must be 12)\n"); }
    if ((win < 1) || (win > 2)) { Rcpp::stop ("Boundary must be either 2 (upper) or 1 (lower)\n"); }

    NumericVector out(length, 0.0);  // Should default to 0s when creating NumericVector, but just in case..

    out = density_PCRM (rts, params, win, step_width);

    return out;
}


// [[Rcpp::export]]
NumericVector dd_IRM (NumericVector rts, NumericVector xj, NumericVector params, int win=1, int method=1)
{
    int length;
    if (rts.length()>1) {
        length = rts.length();
        if ((xj.length() != length) && (xj.length()!=1))  { Rcpp::stop ("rts and xj must have same length or one must be of length 1!\n"); }
        if (xj.length()==1) {
            xj = NumericVector(length, xj[0]);
        }
    } else if (xj.length()>1) {
        length = xj.length();
        if ((rts.length() != length) && (rts.length()!=1))  { Rcpp::stop ("rts and xj must have same length or one must be of length 1!\n"); }
        if (rts.length()==1) {
            rts = NumericVector(length, rts[0]);
        }
    } else {
        length = 1;
    }

    if ((win < 1) || (win > 2)) { Rcpp::stop ("Boundary must be either 2 (upper) or 1 (lower)\n"); }

    NumericVector out(length, 0.0);  // Should default to 0s when creating NumericVector, but just in case..

    double muw = params[win-1];
    double mul = params[2-win];
    double a = params[1+win];
    double b = params[4-win];
    double sigma = params[4];

    if (a==b) {
        double fac = - a/(2*M_PI*sigma*sigma);

        for (int i = 0; i < length; i++) {
            double t = rts[i];
            double sig2t = 2*sigma*sigma*t;
            out[i] = fac /(t*t) *
                exp(-(muw*t+a)*(muw*t+a)/sig2t) *
                (exp(-(xj[i]-(mul*t+a))*(xj[i]-(mul*t+a))/sig2t) -
                exp(-(2*b*mul)/(sigma*sigma))*exp(-(xj[i]-(mul*t-a))*(xj[i]-(mul*t-a))/sig2t));
            if (i % 200 ==0 ) Rcpp::checkUserInterrupt();
        }
    } else {
        if (method==1) {
            double fac = 1/(4*sigma*sigma*M_PI);


            NumericVector C = NumericVector::create(1, -1, -1, 1);
            NumericVector expC1 = NumericVector::create(a, 0, a);
            NumericVector expC2 = NumericVector::create(0, b, b);

            NumericVector expC = -2/(sigma*sigma) * (muw*expC1 + mul*expC2);
            //Rcout << "expC is" << expC << std::endl;
            expC.push_front(0);

            NumericVector Xis = NumericVector::create(a, -a, a, -a,
                                 b, b, -b, -b);
            Xis.attr("dim") = Dimension(4,2);

            for (int i = 0; i < length; i++) {
                double t = rts[i];
                double sig2t = 2*sigma*sigma*t;

                double temp = 0;
                for (int j = 0; j < 4; j++) {
                    double x1tilde = -Xis(j,0)-muw*t;
                    double x2tilde = xj[i] - Xis(j,1) - mul*t;
                    temp = temp + C[j]* x1tilde *
                        exp(expC[j] - (x1tilde)*(x1tilde)/(sig2t) - (x2tilde)*(x2tilde)/(sig2t));
                }
                out[i] = temp*fac/(t*t);
                if (i % 200 ==0 ) Rcpp::checkUserInterrupt();
            }
        } else if (method == 2){
            muw = params[0];
            mul = params[1];
            a = params[2];
            b = params[3];
            double fac = -sigma*sigma/2;
            double muu = muw + mul;
            double muv = muw-mul;

            NumericVector C = NumericVector::create(1, -1, -1, 1);
            NumericVector expC1 = NumericVector::create(0, a, 0, a);
            NumericVector expC2 = NumericVector::create(0, 0, b, b);

            C = C * exp(-2/(sigma*sigma) * (muw*expC1 + mul*expC2));
            //Rcout << "C is" << C << std::endl;

            NumericVector Xis = NumericVector::create(a, -a, a, -a,
                                 b, b, -b, -b);
            Xis.attr("dim") = Dimension(4,2);
            double Sigma = sqrt(2.)*sigma;

            for (int i = 0; i < length; i++) {
                double t = rts[i];
                double fact = -1/(M_PI*Sigma*Sigma*t*t);

                double temp = 0;
                double temp1 = 0;
                for (int j = 0; j < 4; j++) {
                    temp1 = C[j]*((xj[i]-(Xis(j,0)+Xis(j,1))-muu*t)/(Sigma*Sigma) + (3-2*win)*((2*win-3)*xj[i]-(Xis(j,0)-Xis(j,1))-muv*t)/(Sigma*Sigma));
                    temp1 = temp1*exp(-1/(2*t) * ((xj[i]-(Xis(j,0)+Xis(j,1))-muu*t)*(xj[i]-(Xis(j,0)+Xis(j,1))-muu*t)/(Sigma*Sigma) +
                        ((2*win-3)*xj[i]-(Xis(j,0)-Xis(j,1))-muv*t)*((2*win-3)*xj[i]-(Xis(j,0)-Xis(j,1))-muv*t)/(Sigma*Sigma)));
                    temp = temp + temp1;
                }
                out[i] = temp*fact*fac;
                if (i % 200 ==0 ) Rcpp::checkUserInterrupt();
            }
        }
    }

    return out;
}





// [[Rcpp::export]]
NumericVector dd_PCRM (NumericVector rts, NumericVector xj, NumericVector params, int win=1)
{
    int length;
    if (rts.length()>1) {
        length = rts.length();
        if ((xj.length() != length) && (xj.length()!=1))  { Rcpp::stop ("rts and xj must have same length or one must be of length 1!\n"); }
        if (xj.length()==1) {
            xj = NumericVector(length, xj[0]);
        }
    } else if (xj.length()>1) {
        length = xj.length();
        if ((rts.length() != length) && (rts.length()!=1))  { Rcpp::stop ("rts and xj must have same length or one must be of length 1!\n"); }
        if (rts.length()==1) {
            rts = NumericVector(length, rts[0]);
        }
    } else {
        length = 1;
    }



    if ((win < 1) || (win > 2)) { Rcpp::stop ("Boundary must be either 2 (upper) or 1 (lower)\n"); }

    NumericVector out(length, 0.0);  // Should default to 0s when creating NumericVector, but just in case..

    double muw = params[win-1];
    double mul = params[2-win];
    double a = params[1+win];
    double b = params[4-win];
    double sigma = params[4];

    double fac = 1/(3*sqrt(3.)*M_PI*sigma*sigma);

    NumericVector C = NumericVector::create(1, -1, -1, 1, 1, -1);
    NumericVector expC1 = NumericVector::create(0, a, 0, a, a+b, a+b);
    NumericVector expC2 = NumericVector::create(0, 0, b, a+b, b, a+b);

    NumericVector expC = -2/(sigma*sigma) * (muw*expC1 + mul*expC2);

    NumericVector Xis = NumericVector::create(a, -a, a+b, b, -a-b, -b,
                         b, a+b, -b, -a-b, a, -a);
    Xis.attr("dim") = Dimension(6,2);


    for (int i = 0; i < length; i++) {
        double t = rts[i];
        double sig2t = 2*sigma*sigma*t;

        double temp = 0;
        for (int j = 0; j < 6; j++) {
            double x1tilde = -Xis(j,0)-muw*t;
            double x2tilde = xj[i] - Xis(j,1) - mul*t;
            temp = temp + C[j]* (2*x1tilde + x2tilde) *
                exp(expC[j] - (x1tilde)*(x1tilde)/(sig2t) - (x2tilde+0.5*x1tilde)*(x2tilde+0.5*x1tilde)/(0.75*sig2t));
        }
        out[i] = temp*fac/(t*t);
        if (i % 200 ==0 ) Rcpp::checkUserInterrupt();
    }

    return out;
}


// [[Rcpp::export]]
NumericVector r_RM (int n, NumericVector params, double rho, double delta=0.01, double maxT=9)
{
  double sdu, sdv, driftnoise1, driftnoise2; //muu, muv,

  // muu = (params[0]+params[1])*delta;
  // muv = (params[0]-params[1])*delta;

  sdu = sqrt(2*(1+rho)*delta);
  sdv = sqrt(2*(1-rho)*delta);

  double x01, x02, u,v, t, xl;
  int win;

  NumericMatrix out(n, 3);
  for (int i=0; i < n; i++) {
    if (params[6] !=0 ) {
      driftnoise1 = R::rnorm(0, params[6]);
    } else {
      driftnoise1 = 0;
    }
    if (params[7] !=0 ) {
      driftnoise2 = R::rnorm(0, params[7]);
    } else {
      driftnoise2 = 0;
    }
    x01 = params[2]+R::runif(0, params[8]);
    x02 = params[3]+R::runif(0, params[9]);
    t = 0;
    while ((x01 < 0) && (x02 < 0) && (t < maxT)) {
      u = R::rnorm(0, sdu);
      v = R::rnorm(0, sdv);
      x01 = x01 + delta*(driftnoise1 + params[0]) + 0.5*params[4]*(u+v);
      x02 = x02 + delta*(driftnoise2 + params[1]) + 0.5*params[5]*(u-v);
      t += delta;
    }
    if (x01 > 0) {
      if (x02 < x01) {
        win = 1;
        if (x02 < 0) {
          xl = x02;
        } else {
          xl = -1e-24;
        }
      } else {
        win = 2;
        xl = -1e-24;
      }
    } else {
      if (x02 > 0) {
        win = 2;
        xl = x01;
      } else {
        win = 0;
        xl = std::min(x01, x02);
      }
    }
    out( i , 0 ) = t ;
    out( i , 1 ) = win;
    out( i , 2 ) = xl;
    if (i % 200 == 0 ) Rcpp::checkUserInterrupt();

  }
  return out;
}


// [[Rcpp::export]]
NumericMatrix r_WEV (int n, NumericVector params,
                     double delta=0.01, double maxT=9,
                     bool stop_on_error=true)
{
  if (params.length()<13) { Rcpp::stop("Not enough parameters supplied.\n"); }
  NumericMatrix out(n, 6);
  out = RNG_WEV(n,  params, delta, maxT, stop_on_error);

  return out;
}



// [[Rcpp::export]]
NumericVector r_RM_Kiani (int n, NumericVector params, double rho, double Bl, double delta=0.01, double maxT=9)
{
    double a = params[2]; //
    double b = params[3]; //
    double muu = (params[0])*delta;
    double muv = (params[1])*delta;
    double r11 = sqrt((1+sqrt(1-rho*rho))/2)*sqrt(delta)*params[4];
    double r12 = params[4]*params[4]*delta*rho/(2*r11);

    double x01, x02, u,v, t, xl;
    int win;

    NumericMatrix out(n, 3);
    for (int i=0; i < n; i++) {
        x01 = 0;
        x02 = 0;
        t = 0;
        while ((x01 < a) && (x02 < b) && (t < maxT)) {
            u = R::rnorm(0, 1);
            v = R::rnorm(0, 1);
            x01 = x01 + muu + r11*u + r12*v;
            x02 = x02 + muv + r12*u + r11*v;
            if ( x01 < Bl) {
                x01 = Bl - 0.5*(muu + r11*u + r12*v);
            }
            if (x02 < Bl) {
                x02 = Bl - 0.5 * (muv + r12*u + r11*v);
            }
            t += delta;
        }
        if (x01 > a) {
            if (x02 < x01) {
                win = 1;
                if (x02 < b) {
                    xl = x02;
                } else {
                    xl = b-1e-24;
                }
            } else {
                win = 2;
                xl = a-1e-24;
            }
        } else {
            if (x02 > b) {
                win = 2;
                xl = x01;
            } else {
                win = 0;
                xl = -1e-24;
            }
        }

        out( i , 0 ) = t;
        out( i , 1 ) = win;
        out( i , 2 ) = xl;
        if (i % 200 ==0 ) Rcpp::checkUserInterrupt();

    }
    return out;
}


// [[Rcpp::export]]
NumericVector r_LCA (int n, NumericVector params, double delta=0.01, double maxT=9.0)
{
// params: mu1, mu2,
// sig (intensity independent noise component),
// pi (intensity dependency of noise),
// th (threshold), k (leak), beta (inhibition),
// SPV (start point variability), tau (post-dec. accumulation time)
    double mu1 = (params[0]);
    double mu2 = (params[1]);
    double th1 = params[2];
    double th2 = params[3];
    double alpha = params[4]-1;
    double beta = params[5];
    double SPV = params[6];
    double tau = params[7];
    double pi = params[13];
    double sig = params[14];
    double x1, x2, dx1, dx2, t, xl;
    int win;
    double sig1 = sqrt(delta)*sqrt(sig*sig + pi*pi*mu1*mu1);
    double sig2 = sqrt(delta)*sqrt(sig*sig + pi*pi*mu2*mu2);
    NumericMatrix out(n, 5);
    for (int i=0; i < n; i++) {
        x1 = R::runif(0, SPV);
        x2 = R::runif(0, SPV);
        t = 0;
        while ((x1 < th1) && (x2 < th2) && (t < maxT)) {
            dx1 = delta*alpha*x1 - delta*beta*x2 + R::rnorm(mu1*delta, sig1);
            dx2 = delta*alpha*x2 - delta*beta*x1 + R::rnorm(mu2*delta, sig2);
            x1 = std::max(0.0, x1 + dx1);
            x2 = std::max(0.0, x2 + dx2);
            t += delta;
        }
        if (x1 > th1) {
            if (x2 < x1) {
                win = 1;
                if (x2 < th2) {
                    xl = x2;
                } else {
                    xl = th2;
                }
            } else {
                win = 2;
                xl = th1;
            }
        } else {
            if (x2 > th2) {
                win = 2;
                xl = x1;
            } else {
                win = 0;
                xl = std::min(x1, x2);
            }
        }
        if ((win != 0) && (tau>0)) {
            int post_steps = (int)ceil(tau/delta);
            for (int post=1; post<=post_steps; post++) {
                dx1 = delta*alpha*x1 - delta*beta*x2 + R::rnorm(mu1, sig1);
                dx2 = delta*alpha*x2 - delta*beta*x1 + R::rnorm(mu2, sig2);
                x1 = std::max(0.0, x1 + dx1);
                x2 = std::max(0.0, x2 + dx2);
            }
        }
        out( i , 0 ) = t;
        out( i , 1 ) = win;
        out( i , 2 ) = xl;
        out( i , 3 ) = x1;
        out( i , 4 ) = x2;


        if (i % 200 == 0 ) Rcpp::checkUserInterrupt();

    }
    return out;
}


//
// // [[Rcpp::export]]
// NumericVector r_RM2 (int n, NumericVector params, double rho, double delta=0.01, double maxT=9)
// {
//   double sdu, sdv;
//   double mu1 = params[0]*delta;
//   double mu2 = params[1]*delta;
//
//   sdu = sqrt(2*(1+rho)*delta);
//   sdv = sqrt(2*(1-rho)*delta);
//
//   double x01, x02, u1,u2, v, t, xl;
//   int win;
//   double sigrho = 1.0 ;
//   if (rho < 0 ) {
//     rho = -rho;
//     sigrho = -1.0;
//   }
//   NumericMatrix out(n, 3);
//   for (int i=0; i < n; i++) {
//     x01 = params[2];
//     x02 = params[3];
//     t = 0;
//     while ((x01 < 0) && (x02 < 0) && (t < maxT)) {
//       u1 = R::rnorm(0, 1);
//       u2 = R::rnorm(0, 1);
//       v = R::rnorm(0, 1);
//       x01 = x01+mu1+ params[4]*sqrt(delta)*(sqrt((1-rho))*u1+sqrt(rho)*v);
//       x02 = x02+mu2+ params[5]*sqrt(delta)*(sqrt((1-rho))*u2+sqrt(rho)*sigrho*v);
//       t += delta;
//     }
//     if (x01 > 0) {
//       if (x02 < x01) {
//         win = 1;
//         if (x02 < 0) {
//           xl = x02;
//         } else {
//           xl = -1e-24;
//         }
//       } else {
//         win = 2;
//         xl = -1e-24;
//       }
//     } else {
//       if (x02 > 0) {
//         win = 2;
//         xl = x01;
//       } else {
//         win = 0;
//         xl = std::min(x01, x02);
//       }
//     }
//     out( i , 0 ) = t ;
//     out( i , 1 ) = win;
//     out( i , 2 ) = xl;
//     if (i % 200 == 0 ) Rcpp::checkUserInterrupt();
//
//   }
//   return out;
// }
//



// [[Rcpp::export]]
NumericVector r_DDConf (int n, NumericVector params, double delta=0.01, double maxT=9, bool stop_on_error=true)
{
  double a   = params[0];
  double v   = params[1];
  double t0  = params[2];
  double d   = params[3];
  double szr = params[4];
  double sv  = params[5];
  double st0 = params[6];
  double zr  = params[7];
  NumericMatrix out(n, 3);
  double mu, x0, t, decisiont; //conf
  int resp;

  for (int i=0; i < n; i++) {
    mu = R::rnorm(v, sv);
    x0 = a* R::runif(zr-szr/2, zr+szr/2);
    t = 0;
    while ((x0 > 0) && (x0 < a) && (t < maxT)) {
      x0 = x0 + R::rnorm(delta*mu, sqrt(delta));
      t = t + delta;
    }
    if (x0 >= a) {
      resp = 1;
      t = std::max(0.0, t - d/2);
    } else {
      if (x0 <= 0) {
        resp = -1;
        t = std::max(0.0, t + d/2);
      } else {
        resp = 0;
      }
    }
    //conf = 1/sqrt(t);
    decisiont = t;
    t = t  + R::runif(t0-st0/2, t0+st0/2);
    out( i , 0 ) = t;
    out( i , 1 ) = resp;
    out( i , 2 ) = decisiont;
    if (i % 200 ==0 ) Rcpp::checkUserInterrupt();

  }

  return out;

}




// Simulation for the multiple threshold log-normal race model
// [[Rcpp::export]]
NumericVector r_MT_LNR (int n, NumericVector params,
                        NumericVector thresholds)
{
  double mu_d1 = params[0];
  double mu_d2 = params[1];
  double mu_v1 = params[2];
  double mu_v2 = params[3];
  double s_d1  = params[4];
  double s_d2  = params[5];
  double s_v1  = params[6];
  double s_v2  = params[7];
  double rho_d = params[8];
  double rho_v = params[9];
  double t0    = params[10];
  double st0   = params[11];

  int N_thresh = thresholds.length()/2;


  double s_1 = sqrt(s_d1*s_d1+s_v1*s_v1);
  double s_2 = sqrt(s_d2*s_d2+s_v2*s_v2);
  double rho12 = (s_d1*s_d2*rho_d + s_v1*s_v2*rho_v) / (s_1*s_2);


  NumericMatrix out(n, 5);

  double rand_1, rand_2, logT1, logT2, T_dec, logT_ratio;
  int resp, rating;

  for (int i=0; i < n; i++) {
    rand_1 = R::rnorm(0,1);
    rand_2 = R::rnorm(0,1);
    logT1 = mu_d1 - mu_v1 + s_1 * rand_1;
    logT2 = mu_d2 - mu_v2 + s_2 * (rho12 * rand_1 + sqrt(1-rho12*rho12) * rand_2);

    if (logT1 < logT2) {
      resp = 1;
      T_dec = exp(logT1);
      logT_ratio = logT2-logT1;
      rating = 1;
      for (int k=0; k < N_thresh; k++) {
        if (logT_ratio > thresholds[k]) {
          rating = rating + 1;
        }
      }
    } else {
      resp = 2;
      T_dec = exp(logT2);
      logT_ratio = logT1-logT2;
      rating = 1;
      for (int k=0; k < N_thresh; k++) {
        if (logT_ratio > thresholds[ k+N_thresh]) {
          rating = rating + 1;
        }
      }
    }
    out( i , 0 ) = T_dec  + t0 + R::runif(0, st0);   // RT
    out( i , 1 ) = resp;        // response
    out( i , 2 ) = T_dec;       // decision time
    out( i , 3 ) = logT_ratio;  // confidence variable
    out( i , 4 ) = rating;        // discrete confidence rating

    if (i % 200 ==0 ) Rcpp::checkUserInterrupt();

  }

  return out;

}






// // Simulation for the multiple threshold log-normal race model - without any simplifications
// // [[Rcpp::export]]
// NumericVector r_MT_LNR_raw (int n, NumericVector params,
//                         NumericVector thresholds)
// {
//   double mu_d1 = params[0];
//   double mu_d2 = params[1];
//   double mu_v1 = params[2];
//   double mu_v2 = params[3];
//   double s_d1  = params[4];
//   double s_d2  = params[5];
//   double s_v1  = params[6];
//   double s_v2  = params[7];
//   double rho_d = params[8];
//   double rho_v = params[9];
//   double t0    = params[10];
//   double st0   = params[11];
//
//   int N_thresh = thresholds.length()/2;
//
//   NumericMatrix out(n, 5);
//
//   double rand_d1, rand_d2, rand_v1, rand_v2, logT1, logT2, T_dec, logT_ratio;
//   int resp, rating;
//
//   for (int i=0; i < n; i++) {
//     rand_d1 = R::rnorm(0,1);
//     rand_d2 = mu_d2 + s_d2 * (rho_d* rand_d1 + sqrt(1-rho_d*rho_d) * R::rnorm(0,1));
//     rand_d1 = mu_d1 + s_d1 * rand_d1;
//
//     rand_v1 = R::rnorm(0,1);
//     rand_v2 = mu_v2 + s_v2 * (rho_v* rand_v1 + sqrt(1-rho_v*rho_v) * R::rnorm(0,1));
//     rand_v1 = mu_v1 + s_v1 * rand_v1;
//
//     logT1 = rand_d1 - rand_v1;
//     logT2 = rand_d2 - rand_v2;
//
//     if (logT1 < logT2) {
//       resp = 1;
//       T_dec = exp(logT1);
//       logT_ratio = logT2-logT1;
//       rating = 1;
//       for (int k=0; k < N_thresh; k++) {
//         if (logT_ratio > thresholds[k]) {
//           rating = rating + 1;
//         }
//       }
//     } else {
//       resp = 2;
//       T_dec = exp(logT2);
//       logT_ratio = logT1-logT2;
//       rating = 1;
//       for (int k=0; k < N_thresh; k++) {
//         if (logT_ratio > thresholds[ k+N_thresh]) {
//           rating = rating + 1;
//         }
//       }
//     }
//     out( i , 0 ) = T_dec  + t0 + R::runif(0, st0);   // RT
//     out( i , 1 ) = resp;        // response
//     out( i , 2 ) = T_dec;       // decision time
//     out( i , 3 ) = logT_ratio;  // confidence variable
//     out( i , 4 ) = rating;        // discrete confidence rating
//
//
//     if (i % 200 ==0 ) Rcpp::checkUserInterrupt();
//
//   }
//
//   return out;
//
// }


// Probability distribution function for the multiple threshold log-normal race model
// [[Rcpp::export]]
NumericVector d_MT_LNR (NumericVector rts, NumericVector rating,
                        NumericVector params,
                        NumericVector thresholds,
                        int win = 1,
                        double step_width = 0.0001)
{
  if (params.length()!= 12) {Rcpp::stop ("Wrong number of parameters given. (Must be 12)\n"); }
  if ((win < 1) || (win > 2)) { Rcpp::stop ("Boundary must be either 2 or 1!\n"); }

  double mu_d1 = params[0];
  double mu_d2 = params[1];
  double mu_v1 = params[2];
  double mu_v2 = params[3];
  double s_d1  = params[4];
  double s_d2  = params[5];
  double s_v1  = params[6];
  double s_v2  = params[7];
  double rho_d = params[8];
  double rho_v = params[9];
  double t0    = params[10];
  double st0   = params[11];
  if (st0 < EPSILON) {
    st0 = 0;
  }

  // int N_thresh = thresholds.length()/2;
  // NumericVector thresholds_1 = head(thresholds, N_thresh); // first N_thresh entries
  // NumericVector thresholds_2 = tail(thresholds, N_thresh); // last N_thresh entries
  double mu_win, s_win, mu_los, s_los;
  if (win == 1) {
    mu_win = mu_d1 - mu_v1;
    s_win  = sqrt(s_d1*s_d1 + s_v1*s_v1);
    mu_los = mu_d2 - mu_v2;
    s_los  = sqrt(s_d2*s_d2 + s_v2*s_v2);
  } else {
    mu_win = mu_d2 - mu_v2;
    s_win  = sqrt(s_d2*s_d2 + s_v2*s_v2);
    mu_los = mu_d1 - mu_v1;
    s_los  = sqrt(s_d1*s_d1 + s_v1*s_v1);
  }
  double rho_12 = (s_d1*s_d2*rho_d + s_v1*s_v2*rho_v) / (s_win*s_los);
  double th1, th2;
  int length = rts.length();

  NumericVector out(length, 0.0);  // Should default to 0s when creating NumericVector, but just in case..
  for (int i=0; i < length; i++) {
    if (rts[i]-t0 <= 0 ) {
      out[i] = 0;
    } else {
      th2= thresholds[rating[i]];
      th1= thresholds[rating[i]-1];
      out[i] = density_MTLNR (rts[i]-t0, th1, th2,
                              mu_win, mu_los,
                              s_win, s_los, rho_12,
                              t0, st0, step_width);
    }
    if (i % 200 ==0 ) Rcpp::checkUserInterrupt();
  }
  return out;

}



// End CPP
