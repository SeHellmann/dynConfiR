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

  NumericMatrix out(n, 6);

  double a   =    params[0];
  double v   =    params[1];
  double t0  =    params[2];
  double d   =    params[3];
  double szr =    params[4];
  double sv  =    params[5];
  double st0 =    params[6];
  double zr  =    params[7];
  double tau =    params[8];
  double lambda =  params[11];
  double w   =    params[12];
  double muvis =  params[13];
  double sigvis = params[14];
  double svis =   params[15];

  double mu, x0, t, conf, vis;
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
    } else {
      if (x0 <= 0) {
        resp = -1;
      } else {
        resp = 0;
      }
    }
    //if (w == 1) {
    if (tau > 0) {
      conf = resp*(x0 + R::rnorm(tau*mu, sqrt(tau)) - a*zr);
    } else {
      conf = resp*(x0 - a*zr);
    }
    // save response, response time and state of evidence accumulator
    out( i , 0 ) = std::max(0.0, t - resp*d/2)  + R::runif(t0-st0/2, t0+st0/2);;
    out( i , 1 ) = resp;
    out( i , 3 ) = conf; // evidence term

    vis = R::rnorm((tau+t)*muvis, sqrt(svis*svis*(tau+t)+(t+tau)*(t+tau)*sigvis*sigvis));

    if (lambda >0) {
      conf = (w*conf + (1-w)*vis)/pow(t+tau, lambda);
    } else {
      conf = (w*conf + (1-w)*vis);
    }
    // save final confidence and value of visibility process
    out( i , 2 ) = conf;
    out( i , 4 ) = vis;
    out( i , 5 ) = mu;

    if (i % 200 ==0 ) Rcpp::checkUserInterrupt();

  }

  return out;
}

#endif // RNG_WEV_H
