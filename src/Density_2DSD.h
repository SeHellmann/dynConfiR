/* Density_2DSD.h - Functions for PDF calculation of the 2DSD model
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

#ifndef DENSITY_2DSD_H
#define DENSITY_2DSD_H

using namespace Rcpp;

#define EPSILON 1e-6

// Forward declarations
double g_minus_2DSD (double t, NumericVector params);

static double integral_t0_g_minus_2DSD (double t, NumericVector params);
static double integral_z_g_minus_2DSD  (double t, NumericVector params);
static double integral_v_g_minus_2DSD  (double t, double zr, NumericVector params);

static double g_minus_no_var_2DSD     (double t, double a, double zr, double v, double tau, double th1, double th2, double lambda);
static double g_minus_small_time_2DSD (double t, double zr, int N);
static double g_minus_large_time_2DSD (double t, double zr, int N);

// TODO: Make sure these function names are accurate
static double integrate_z_over_t_2DSD  (NumericVector params, double a, double b, double step_width);
static double integrate_v_over_zr_2DSD (NumericVector params, double a, double b, double t, double step_width);


// Main calls
NumericVector density_2DSD (NumericVector rts, NumericVector params, int boundary, int stopon0)
{
    int length = rts.length();
    NumericVector out(length);
    if (stopon0==1) {
      if (boundary == 1) {
        params[7] = 1- params[7]; // z -> 1 - z
        params[1] = - params[1]; // v  -> - v
        params[3] = - params[3]; // d  -> - d
        for (int i = 0; i < length; i++) {
          out[i] =  g_minus_2DSD(rts[i], params);
          if (out[i]==0) break;
          if (i % 200 ==0 ) Rcpp::checkUserInterrupt();
        }
      } // Calc upper
      else {
        for (int i = 0; i < length; i++) {
          out[i] = -g_minus_2DSD(rts[i], params);
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
          out[i] =  g_minus_2DSD(rts[i], params);
          if (i % 200 ==0 ) Rcpp::checkUserInterrupt();
        }
      } // Calc upper
      else {
        for (int i = 0; i < length; i++) {
          out[i] = -g_minus_2DSD(rts[i], params);
          if (i % 200 ==0 ) Rcpp::checkUserInterrupt();
        }
      } // Calc lower
    }


    return out;
}

double g_minus_2DSD(double t, NumericVector params)
{
    return integral_t0_g_minus_2DSD (t - params[2] - 0.5*params[3], params);
}



static double integral_t0_g_minus_2DSD (double t, NumericVector params)
{
    double res;

    if (params[6] < params[15]) // 170501   was == 0)
    {
        res = integral_z_g_minus_2DSD (t, params);
    }
    else
    {
        res = integrate_z_over_t_2DSD(params,
                        t - .5*params[6],
                        t + .5*params[6], params[12]) / params[6];
    }

    return res;
}


static double integral_z_g_minus_2DSD (double t, NumericVector params)
{
    double res;

    if (t <= 0) return 0;

    if (params[4] < params[14])
    {
        res = integral_v_g_minus_2DSD (t, params[7], params);
    }
    else
    {
        res = integrate_v_over_zr_2DSD(params, params[7] - .5*params[4], params[7] + .5*params[4],
                                  t, params[13]) / params[4];
    }
    return res;
}


static double integral_v_g_minus_2DSD (double t, double zr, NumericVector params)
{
    double a = params[0];
    double v = params[1];
    double sv = params[5];
    double tau = params[8];
    double th2 = params[10];
    double th1 = params[9];
    double lambda = params[11];

    if (sv == 0)
    {
        return g_minus_no_var_2DSD(t, a, zr, v, tau, th1, th2, lambda);
    }

    int N_small, N_large;
    double sv2t, mean_conf, sd_conf,diff_normal, simple, factor, eps;

    double ta = t/(a*a);

    sv2t = sv*sv*t + 1;
    factor = 1 / (a*a * sqrt(sv2t)) * exp(-0.5 * (v*v*t + 2*a*zr*v - a*a * zr*zr*sv*sv  ) / (sv2t)); //sqrt(t+tau)

    mean_conf = - (tau*v - a*zr*(sv*sv*(tau+t)+1)) / (sv2t); // ((v - sv*sv*a*zr)*tau) / sv2t;
    sd_conf = sqrt( tau*(sv*sv*tau+sv2t)/(sv2t));
    diff_normal = 0.5*(erf((th2*pow(t+tau, lambda) - mean_conf) / (M_SQRT2 * sd_conf)) -
                        erf((th1*pow(t+tau, lambda) - mean_conf) / (M_SQRT2 * sd_conf)));

    // factor = exp(-a*zr*v - 0.5*v*v*t) / (a*a); //sqrt(t+tau)
    // diff_normal = 0.5*(erf((th2*pow(t+tau, lambda) + (tau*v - a*zr)) / (M_SQRT2 * sqrt(tau))) -
    //   erf((th1*pow(t+tau, lambda) + (tau*v - a*zr)) / (M_SQRT2 * sqrt(tau))));
    //
    //

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
        simple = g_minus_small_time_2DSD(t/(a*a), zr, N_small);
    }
    else
    {
        simple = g_minus_large_time_2DSD(t/(a*a), zr, N_large);
    }
    return factor * diff_normal * simple;
}


static double g_minus_no_var_2DSD(double t, double a, double zr, double v,
                                   double tau, double th1, double th2, double lambda)
{
    int N_small, N_large;
    double diff_normal, simple, factor, eps;
    double ta = t/(a*a);

    factor = exp(-a*zr*v - 0.5*v*v*t) / (a*a); //sqrt(t+tau)
    diff_normal = 0.5*(erf((th2*pow(t+tau, lambda) + (tau*v - a*zr)) / (M_SQRT2 * sqrt(tau))) -
      erf((th1*pow(t+tau, lambda) + (tau*v - a*zr)) / (M_SQRT2 * sqrt(tau))));

    /*
    if (std::isinf(factor))
    {
      return 0;
    }
     */
  /*  if (diff_normal == 0)
    {
      return 0;
    }    */

    eps = EPSILON / factor;

    N_large = (int)ceil (1/ (M_PI*sqrt(t)));
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
        simple = g_minus_small_time_2DSD(t/(a*a), zr, 2*N_small);
    }
    else
    {
        simple = g_minus_large_time_2DSD(t/(a*a), zr, 2*N_large);
    }
    return factor * diff_normal * simple;
}


static double g_minus_small_time_2DSD(double t, double zr, int N)
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

static double g_minus_large_time_2DSD(double t, double zr, int N)
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
static double integrate_z_over_t_2DSD (NumericVector params, double a, double b, double step_width)
{
    double width = b-a;
    int N = std::max(4, (int) (width / step_width));
    double step = width / N; // std::max(width / N, EPSILON);
    double x;
    double result = 0;

    for(x = a+0.5*step; x < b; x += step)
    {
        result += step * integral_z_g_minus_2DSD(x, params);
    }
    return result;
}

static double integrate_v_over_zr_2DSD (NumericVector params, double a, double b, double t, double step_width)
{
    double width = b-a;
    int N = std::max(4, (int) (width / step_width));
    double step = width / N; // std::max(width / N, EPSILON);
    double x;
    double result = 0;

    for(x = a+0.5*step; x < b; x += step)
    {
      result += step * integral_v_g_minus_2DSD (t, x, params);
    }
    return result;
}

#endif // DENSITY_2DSD_H
