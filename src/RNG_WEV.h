/* RNG_WEV.h - Functions for generating random trials in the WEV model
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

#ifndef RNG_WEV_H
#define RNG_WEV_H

using namespace Rcpp;

NumericMatrix RNG_WEV (int n, NumericVector params, double delta=0.01,
                     double maxT=9, bool stop_on_error=true)
{

  NumericMatrix out(n, 3);

  double a   = params[0];
  double v   = params[1];
  double t0  = params[2];
  double d   = params[3];
  double szr = params[4];
  double sv  = params[5];
  double st0 = params[6];
  double zr  = params[7];
  double tau = params[8];
  double w   = params[11];
  double muvis = params[12];
  double sigvis = params[13];
  double svis = params[14];

  double mu, x0, t, conf, evid, vis;
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
    if (w == 1) {
      conf = x0 + R::rnorm(tau*mu, sqrt(tau));
    } else {
      evid = R::rnorm(tau*mu*resp, sqrt(tau));
      vis = R::rnorm((tau+t)*muvis, sqrt(svis*svis*(tau+t)+(t+tau)*(t+tau)*sigvis*sigvis));
      conf = w*evid + (1-w)*vis;
    }
    t = t  + R::runif(t0-st0/2, t0+st0/2);
    out( i , 0 ) = t;
    out( i , 1 ) = resp;
    out( i , 2 ) = conf;
    if (i % 200 ==0 ) Rcpp::checkUserInterrupt();

  }

  return out;
}

#endif // RNG_WEV_H
