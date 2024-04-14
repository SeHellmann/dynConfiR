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

#ifndef VALIDATEPARAMETERS_H
#define VALIDATEPARAMETERS_H

using namespace Rcpp;

bool ValidateParams (NumericVector params, bool print)
{
  bool valid = true;
  // Report error, if:
  // a <= 0
  if (params[0] <= 0)                 { valid = false; if (print) Rcpp::Rcout << "error: invalid parameter a = " << params[0] << std::endl;  }
  // z < 0 OR z > 1
  if (params[7] < 0 || params[7] > 1) { valid = false; if (print) Rcpp::Rcout << "error: invalid parameter zr = " << params[7] << std::endl; }
  // st0 < 0
  if (params[6] < 0)                  { valid = false; if (print) Rcpp::Rcout << "error: invalid parameter st0 = " << params[6] << std::endl; }
  // szr < 0
  if (params[4] < 0)                  { valid = false; if (print) Rcpp::Rcout << "error: invalid parameter szr = " << params[4] << std::endl; }
  // sv < 0
  if (params[5] < 0)                  { valid = false; if (print) Rcpp::Rcout << "error: invalid parameter sv = " << params[5] << std::endl; }
  // t0 < 0.5*(d+st0)
  if (params[2] - fabs(0.5*params[3]) - 0.5*params[6] < 0) { valid = false; if (print) Rcpp::Rcout << "error: invalid parameter combination t0 = " << params[2] << ", d = " << params[3] << ", st0 =" << params[6] << std::endl; }
  // z - sz/2 < 0 OR z + sz/2 > 1
  if (params[7] - 0.5*params[4] < 0) { valid = false; if (print) Rcpp::Rcout << "error: invalid parameter combination zr = " << params[7] << ", szr = " << params[4] << std::endl;}
  if (params[7] + 0.5*params[4] > 1) { valid = false; if (print) Rcpp::Rcout << "error: invalid parameter combination zr = " << params[7] << ", szr = " << params[4] << std::endl;}

  if (params.size()<= 10) { // DDConf model
    // th2 < th1
    if (params[9] < params[8])      { valid = false; if (print) Rcpp::Rcout << "error: invalid parameter combination th1 = " << params[8] << ", th2 = " << params[9] << std::endl;}
  } else { // 2DSDT and dynaViTE
    // tau < 0
    if (params[8] < 0)              { valid = false; if (print) Rcpp::Rcout << "error: invalid parameter tau = " << params[8] << std::endl;}
    // th2 < th1
    if (params[10] < params[9])     { valid = false; if (print) Rcpp::Rcout << "error: invalid parameter combination th1 = " << params[9] << ", th2 = " << params[10] << std::endl;}
    // lambda < 0
    if (params[11] < 0)             { valid = false; if (print) Rcpp::Rcout << "error: invalid parameter lambda = " << params[11] <<  std::endl; }
    if (params.size() > 12) { // only dynaViTE
      // w < 0 OR w > 1
      if (params[12]<0 || params[12] > 1)             { valid = false; if (print) Rcpp::Rcout << "error: invalid parameter w = " << params[12] << ", allowed: w in [0,1]" <<  std::endl; }
      // sigvis < 0
      if (params[14] < 0)               { valid = false; if (print) Rcpp::Rcout << "error: invalid parameter sigvis = " << params[14] <<  std::endl; }
      // svis <= 0
      if (params[15] <= 0)                { valid = false; if (print) Rcpp::Rcout << "error: invalid parameter svis = " << params[15] <<  std::endl; }
    }
  }
  return valid;
}

#endif // VALIDATEPARAMETERS_H
