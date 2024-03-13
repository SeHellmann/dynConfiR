/* Density_WEVmu.h - Functions for PDF calculation in dynWEV model
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

#ifndef DENSITY_WEVmu_H
#define DENSITY_WEVmu_H

using namespace Rcpp;

#define EPSILON 1e-6

// Forward declarations

double g_minus_WEVmu (double t, NumericVector params);

static double integral_t0_g_minus_WEVmu (double t, NumericVector params);
static double integral_z_g_minus_WEVmu  (double t, NumericVector params);
static double integral_v_g_minus_WEVmu  (double t, double zr, NumericVector params);

static double g_minus_small_time_WEVmu (double t, double zr, int N);
static double g_minus_large_time_WEVmu (double t, double zr, int N);

// TODO: Make sure these function names are accurate
static double integrate_z_over_t_WEVmu  (NumericVector params, double a, double b, double step_width);
static double integrate_v_over_zr_WEVmu (NumericVector params, double a, double b, double t, double step_width);


// Main calls
NumericVector density_WEVmu (NumericVector rts, NumericVector params, int boundary, int stopon0)
{
    int length = rts.length();
    NumericVector out(length);
    if (stopon0 == 1) {
      if (boundary == 1) {
        params[7] = 1- params[7]; // z -> 1 - z
        params[1] = - params[1]; // v  -> - v
        params[3] = - params[3]; // d  -> - d
        for (int i = 0; i < length; i++) {
          out[i] =  g_minus_WEVmu(rts[i], params);
          if (out[i]==0) break;
          if (i % 200 ==0 ) Rcpp::checkUserInterrupt();
        }
      } // Calc upper
      else {
        for (int i = 0; i < length; i++) {
          out[i] = -g_minus_WEVmu(rts[i], params);
          if (out[i]==0) break;
          if (i % 200 ==0 ) Rcpp::checkUserInterrupt();
        }
      } // Calc lower
    } else {
      if (boundary == 1) {
        params[7] = 1- params[7]; // z -> 1 - z
        params[1] = - params[1]; // v  -> - v
        params[3] = - params[3]; // d  -> - d
        for (int i = 0; i < length; i++) {
          out[i] =  g_minus_WEVmu(rts[i], params);
          if (i % 200 ==0 ) Rcpp::checkUserInterrupt();
          }
        } // Calc upper
      else {
        for (int i = 0; i < length; i++) {
          out[i] = -g_minus_WEVmu(rts[i], params);
          if (i % 200 ==0 ) Rcpp::checkUserInterrupt();
        }
      } // Calc lower
    }
    return out;
}

double density_WEVmu2 (double rt, NumericVector params, int boundary, int stopon0)
{
  if (boundary == 1) {
      params[7] = 1- params[7]; // z -> 1 - z
      params[1] = - params[1]; // v  -> - v
      params[3] = - params[3]; // d  -> - d
    } // Calc upper
  return -g_minus_WEVmu(rt, params);  // Calc lower
}


double g_minus_WEVmu(double t, NumericVector params)
{
    return integral_t0_g_minus_WEVmu (t - params[2] - 0.5*params[3], params);
}



static double integral_t0_g_minus_WEVmu (double t, NumericVector params)
{
    double res;

    if (params[6] < params[19]) // 170501   was == 0)
    {
        res = integral_z_g_minus_WEVmu (t, params);
    }
    else
    {
        res = integrate_z_over_t_WEVmu(params,
                                       t - .5*params[6],
                                       t + .5*params[6], params[16]) / params[6];
    }

    return res;
}


static double integral_z_g_minus_WEVmu (double t, NumericVector params)
{
    double res;

    if (t <= 0) return 0;

    if (params[4] < params[18])
    {
        res = integral_v_g_minus_WEVmu (t, params[7], params);
    }
    else
    {
        res = integrate_v_over_zr_WEVmu(params, params[7] - .5*params[4], params[7] + .5*params[4],
                                        t, params[17]) / params[4];
    }
    return res;
}


static double integral_v_g_minus_WEVmu (double t, double zr, NumericVector params)
{
  double a = params[0];
  double v = params[1];
  double sv = params[5];
  double tau = params[8];
  double th1 = params[9];
  double th2 = params[10];
  double lambda = params[11];
  double w = params[12];
  double muvis = params[13];
  double sigvis2 = params[14]*params[14];
  double svis2 = params[15]*params[15];



    // Compute the factor in front of integral and series
    int N_small, N_large;
    double sv2t, simple, factor, eps, ttau, Mu, Sigma, int_c;


    double ta = t/(a*a);

    sv2t = sv*sv*t + 1;
    factor = 1 / (a*a * sqrt(sv2t)) * exp(-0.5 * (v*v*t + 2*a*zr*v - a*a * zr*zr*sv*sv  ) / (sv2t));

    //if (std::isinf(factor)||(factor==0)) { return 0; }
    if (factor==0) { return 0; }

    // Compute the integral w.r.t. the confidence variable c:
    ttau = t+tau;
    Mu = (ttau*(1-w)*muvis - w*(tau*v-a*zr*(sv*sv*ttau+1))/(sv2t)) ;
    Sigma = sqrt(( w*w*tau* (1+ tau*sv*sv/sv2t) + (1-w)*(1-w)*(ttau*svis2 + ttau*ttau*sigvis2)));
    if (lambda > 0)
    {
      int_c = Phi((th2*pow(ttau, lambda)-Mu)/Sigma) - Phi((th1*pow(ttau, lambda)-Mu)/Sigma);
    }
    else
    {
      int_c = Phi((th2-Mu)/Sigma) - Phi((th1-Mu)/Sigma);
    }

    // Approximate the series in the function f:
    eps = EPSILON / factor;


    N_large = (int)ceil(1 / (M_PI*sqrt(t)));
    if (M_PI*ta*eps < 1)
    {
        N_large = std::max(N_large, (int)ceil(sqrt(-2*log(M_PI*ta*eps) / (M_PI*M_PI*ta))));
    }

    if (2*sqrt(2*M_PI*ta)*eps < 1)
    {
        N_small = (int)ceil(fmax(sqrt(ta) + 1, 2 + sqrt(-2*ta*log(2*eps*sqrt(2*M_PI*ta)))));
    }
    else
    {
        N_small = 2;
    }

    if (N_small < N_large)
    {
        simple = g_minus_small_time_WEVmu(t/(a*a), zr, N_small);
    }
    else
    {
        simple = g_minus_large_time_WEVmu(t/(a*a), zr, N_large);
    }
    return factor * int_c * simple;
}


static double g_minus_small_time_WEVmu(double t, double zr, int N)
{
    int i;
    double sum = 0;
    double d;

    for(i = -N/2; i <= N/2; i++)
    {
        d = 2*i + zr;
        sum += exp(-d*d / (2*t)) * d;
    }

    return sum / sqrt(2*M_PI*t*t*t);
}

static double g_minus_large_time_WEVmu(double t, double zr, int N)
{
    int i;
    double sum = 0;
    double d;

    for(i = 1; i <= N; i++)
    {
        d = i * M_PI;
        sum += exp(-0.5 * d*d * t) * sin(d*zr) * i;
    }

    return sum * M_PI;
}

// CONVERSION NOTE: Simplest way to deal with the integrate function is to remove
//                  the clever recursiveness and instead (ugh) duplicate code
static double integrate_z_over_t_WEVmu (NumericVector params, double a, double b, double step_width)
{
    double width = b-a;
    int N = std::max(4, (int) (width / step_width));
    double step = std::max(width / N, EPSILON);
    double x;
    double result = 0;

    for(x = a+0.5*step; x < b; x += step)
    {
        result += step * integral_z_g_minus_WEVmu(x, params);
    }
    return result;
}

static double integrate_v_over_zr_WEVmu (NumericVector params, double a, double b, double t, double step_width)
{
    double width = b-a;
    int N = std::max(4, (int) (width / step_width));
    double step = std::max(width / N, EPSILON);
    double x;
    double result = 0;

    for(x = a+0.5*step; x < b; x += step)
    {
      result += step * integral_v_g_minus_WEVmu (t, x, params);
    }
    return result;
}


#endif // DENSITY_WEVmu_H
