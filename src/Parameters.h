/* Parameters.h - A class to contain the model parameters and precision tuning
 *   (originally in parameters.c and precision.c)
 *
 * Copyright (C) 2006  Jochen Voss, Andreas Voss.
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation; either version 2 of the
 * License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 */

#ifndef PARAMETERS_H
#define PARAMETERS_H

// Note: Parameters class now includes precision constants

// Indices for packed parameters array
#define PARAM_a    0
#define PARAM_v    1
#define PARAM_t0   2
#define PARAM_d    3
#define PARAM_szr  4
#define PARAM_sv   5
#define PARAM_st0  6
#define PARAM_zr   7
#define PARAM_tau  8
#define PARAM_th1  9
#define PARAM_th2  10
#define PARAM_w    11
#define PARAM_muvis  12
#define PARAM_sigvis 13
#define PARAM_svis   14


class Parameters
{
public:
    double a;     // Boundary separation
    double v;     // Mean of the drift
    double t0;    // Non-decision time
    double d;     // Difference between boundaries of non-decision time
    double szr;   // width of zr distribution
    double sv;    // standard deviation of v distribution
    double st0;   // width of t0 distribution
    double zr;    // Mean of diffusion starting point relative to boundaries
    double tau;   // Postdecisional accumulation time for confidence
    double th1;   // Lower bound for confidence interval
    double th2;   // Upper bound for confidence interval
    double thD;   // State of the process at decision time (a for upper; 0 for lower boundary)
    double w;     // weight on decision evidence; (1-w) is weight on visibility accumulator
    double muvis;   // mean drift in the visibility accumulation process
    double sigvis; // variability in the drift rate of the visibility process
    double svis;   // variability in the visibility accumulation process
    double q_WEV; // frequently used constant in WEV version 2 and 3

    // Precision constants set by SetPrecision()
    double  TUNE_INT_T0;
    double  TUNE_INT_Z;

    double  TUNE_SV_EPSILON; // CONVERSION NOTE: See below in SetPrecision()
    double  TUNE_SZ_EPSILON; // CONVERSION NOTE: See below in SetPrecision()
    double  TUNE_ST0_EPSILON; // CONVERSION NOTE: See below in SetPrecision()

public:
  // Construct the object from the passed in params
    Parameters (NumericVector params, double precision)   // Construct parameters for the 2DSD Model
    {
        a   = params[PARAM_a];
        v   = params[PARAM_v];
        t0  = params[PARAM_t0];
        d   = params[PARAM_d];
        szr = params[PARAM_szr];
        sv  = params[PARAM_sv];
        st0 = params[PARAM_st0];
        zr  = params[PARAM_zr];
        tau = params[PARAM_tau];
        th1 = params[PARAM_th1];
        th2 = params[PARAM_th2];
        thD = 0; // The default value is 0, because we use g_minus as default function (which corresponds to the lower boundary)
        q_WEV = 0;

        //Rcpp::Rcout << "note: length of params vector " << params.size() << std::endl;

        if (params.size()<= 11)
        {

          w = 0;
          muvis = 0;
          sigvis = 1;
          svis = 1;
        }

        if (params.size() > 11)    // This switches between parameters for dWEV model or simpler 2DSD model
        {
          w   = params[PARAM_w];
          muvis = params[PARAM_muvis];
          sigvis=1;
          svis = 1;
          if (params.size()>13)
          {
            sigvis = params[PARAM_sigvis];
            svis = params[PARAM_svis];
            q_WEV = -(1-w)/(w*tau);
            th2 = -th1/(w*tau);
            th1 = -params[PARAM_th2]/(w*tau);
          }
        }

        SetPrecision (precision);
    }


    bool ValidateParams_2DSD (bool print)
    {
        bool valid = true;

        if (a <= 0)                         { valid = false; if (print) Rcpp::Rcout << "error: invalid parameter a = " << a << std::endl;  }
        if (szr < 0 || szr > 1)             { valid = false; if (print) Rcpp::Rcout << "error: invalid parameter szr = " << szr << std::endl; }
        if (st0 < 0)                        { valid = false; if (print) Rcpp::Rcout << "error: invalid parameter st0 = " << st0 << std::endl; }
        if (sv < 0)                         { valid = false; if (print) Rcpp::Rcout << "error: invalid parameter sv = " << sv << std::endl; }
        if (t0 - fabs(0.5*d) - 0.5*st0 < 0) { valid = false; if (print) Rcpp::Rcout << "error: invalid parameter combination t0 = " << t0 << ", d = " << d << ", st0 =" << st0 << std::endl; }
        if (zr - 0.5*szr <= 0)              { valid = false; if (print) Rcpp::Rcout << "error: invalid parameter combination zr = " << zr << ", szr = " << szr << std::endl;}
        if (zr + 0.5*szr >= 1)              { valid = false; if (print) Rcpp::Rcout << "error: invalid parameter combination zr = " << zr << ", szr = " << szr << std::endl;}
        if (tau < 0)                        { valid = false; if (print) Rcpp::Rcout << "error: invalid parameter tau = " << tau << std::endl;}
        if (th2 < th1)                      { valid = false; if (print) Rcpp::Rcout << "error: invalid parameter combination th1 = " << th1 << ", th2 = " << th2 << std::endl;}
        if (w<0 || w >= 1)                  { valid = false; if (print) Rcpp::Rcout << "error: invalid parameter w = " << w << ", allowed: w in [0,1)" <<  std::endl; }
        if (sigvis < 0)                      { valid = false; if (print) Rcpp::Rcout << "error: invalid parameter sigvis = " << sigvis <<  std::endl; }
        if (svis <= 0)                        { valid = false; if (print) Rcpp::Rcout << "error: invalid parameter svis = " << svis <<  std::endl; }

        return valid;
    }

private:
    void  SetPrecision (double p)
    {
        TUNE_INT_T0 = 0.0089045 * exp(-1.037580*p);
        TUNE_INT_Z  = 0.0508061 * exp(-1.022373*p);  // Wie kommt man dazu?? --> Voss

        // CONVERSION NOTE:
        //     These have been added to optimise code paths by treating very small variances as 0
        //     e.g. with precision = 3, sv or sz values < 10^-5 are considered 0
        TUNE_SZ_EPSILON = pow (10, -(p+2.0));     // Used by ddiffusion and pdiffusion
        TUNE_ST0_EPSILON = pow (10, -(p+2.0));     // Used by ddiffusion
    }
};

#endif // PARAMETERS_H
