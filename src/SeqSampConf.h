/* SeqSampConf.h - Main header file for the RCpp implementation of the
 * sequential sampling confidence models contained in dynConfiR
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

#ifndef SeqSampConf_H
#define SeqSampConf_H

#include <Rcpp.h>
using namespace Rcpp;
#include <assert.h>

#define MAX_INPUT_VALUES 1e+6

// Include all .hpp files
// Note: This is bad organisation, but Rcpp (and especially RStudio's "sourceCpp" don't seem to handle
//       projects with multiple .cpp files)


// Used by both PDF and CDF (also includes tuning params) also by 2DSD and WEV implementation
#include "Parameters.h"
#include "ValidateParameters.h"


// While not enforced, this is the global parameters Singleton
//   To be created and freed in the other calls in SeqSampConf.cpp
Parameters *g_Params;

#define BOUNDARY_LOWER 0
#define BOUNDARY_UPPER 1

// static double phi(double x, double m, double s);  //pdf of normal distribution with mean m and sd s
static double Phi(double x); //CDF of standard normal distribution

// PDF for 2DSD
#include "Density_2DSD.h"
// PDF and RNG for WEV version 1
#include "Density_WEVmu.h"
#include "RNG_WEV.h"

// PDF and RNG for IRM
#include "Density_IRM.h"
#include "Density_IRM2.h"
// PDF and RNG for PCRM
#include "Density_PCRM.h"
// PDF and RNG for DDMConf
#include "Density_DDMConf.h"


static double Phi (double x)  // The distribution function of the standard normal distribution.
{
        return  0.5*(1+erf (x/M_SQRT2));
}


//
// static double phi(double x, double m, double s) // Pdf of normal distribution with mean m and sd s
// {
//         static const float inv_sqrt_2pi = 0.3989422804014327;
//         float expon = (x - m) / s;
//
//         return inv_sqrt_2pi / s * std::exp(-0.5 * expon*expon);
// }


#endif // SeqSampConf_H
