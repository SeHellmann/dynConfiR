/* Density_MTLNR.h - Functions for PDF calculation in the multiple-threshold log-normal accumulator model
 *
 * Copyright (C) 2025 Sebastian Hellmann.
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

#ifndef DENSITY_MTLNR_H
#define DENSITY_MTLNR_H

using namespace Rcpp;

#define EPSILON 1e-6

static double densMTLNR_const_t (double rt, double th1, double th2,
                                 double mu_win, double mu_los,
                                 double s_win, double s_los,
                                 double rho_12, double denom_erf) {
    double logT = log(rt);
    double fT = exp(-(logT-mu_win)*(logT-mu_win)/(2*s_win*s_win))/(rt*s_win);
    double num_erf2 = (th2+ logT - (mu_los + s_los/s_win * rho_12*(logT-mu_win)));
    double num_erf1 = (th1+ logT - (mu_los + s_los/s_win * rho_12*(logT-mu_win)));
    // Rcout << "logT : " << logT << "\n";
    // Rcout << "fT : " << fT << "\n";
    // Rcout << "num_erf2 : " << num_erf2 << "\n";
    // Rcout << "num_erf1 : " << num_erf1 << "\n";
    // Rcout << "erf(2) : " << num_erf2 << "\n";

    double res = fT *
        (erf(num_erf2/denom_erf) - erf(num_erf1/denom_erf));
    return res;
}

// Main calls
double density_MTLNR (double rt,
                             double th1, double th2,
                             double mu_win, double mu_los,
                             double s_win, double s_los, double rho_12,
                             double t0, double st0,
                             double step_width = 0.0001)
{
    double res;
    static const float inv_sqrt_2pi = 0.3989422804014327;
    double denom_erf = sqrt(2) * s_los*sqrt(1-rho_12*rho_12);
    //Rcout << "The value of denom_erf : " << denom_erf << "\n";
    if (th1 > th2) {
        Rcout << "Stop, because th1: " << th1 << " was bigger than th2:" << th2 << "\n";
        return 0;
    }

    if (st0==0) {
        //Rcout << "Yey, st0=0" << "\n";
        res = inv_sqrt_2pi * 0.5 * densMTLNR_const_t (rt, th1, th2,
                                                      mu_win, mu_los, s_win, s_los,
                                                      rho_12, denom_erf);
    } else {
        int nsteps = std::max(4, (int) (st0/step_width));
        double dt = st0 / nsteps;
        double result = 0;
        for(double x = rt-st0; x < rt; x += dt)
        {
            if ( x > 0)
            {
                result +=  densMTLNR_const_t ( x, th1, th2,
                                mu_win, mu_los, s_win, s_los,
                                rho_12, denom_erf);
            }
        }
        res = inv_sqrt_2pi * 0.5 * dt * result / st0;
    }
    return res;
}

#endif // DENSITY_MTLNR_H
