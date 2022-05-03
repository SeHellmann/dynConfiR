#Actual function for integrating and calculating the probabilities
preddist <- function(row, thetas_lower, thetas_upper, paramDf,V, SV, model, simult_conf,
                     precision=1e-5,
                     maxrt=15, subdivisions=100L,
                     stop.on.error = FALSE,
                     .progress = TRUE, pb=NULL) {
  if (model != "2DSD") {
    vth1 <- ifelse(row$response =="upper", thetas_upper[row$rating], thetas_lower[(row$rating)])
    vth2 <- ifelse(row$response =="upper", thetas_upper[(row$rating+1)], thetas_lower[(row$rating+1)])
  } else {
    vth1 <- ifelse(row$response =="upper", thetas_upper[row$rating], rev(thetas_lower)[(row$rating+1)])
    vth2 <- ifelse(row$response =="upper", thetas_upper[(row$rating+1)], rev(thetas_lower)[(row$rating)])
  }


  integrand <- with(paramDf, switch(which(model== c("dynWEV", "2DSD")),
                                    function(t) return(dWEV(t, vth1,vth2,
                                                                  response=as.character(row$response), tau=tau, a=a,
                                                                  v = (-1)^(row$stimulus=="lower")*V[row$condition],
                                                                  t0 = t0, z = z, sz = sz, st0=st0,
                                                                  sv = SV[row$condition], w=w, svis=svis, sigvis=sigvis,
                                                            simult_conf = simult_conf,
                                                            z_absolute = FALSE, precision = precision)),
                                    function(t) return(d2DSD(t, vth1,vth2,
                                                             response=as.character(row$response), tau=tau, a=a,
                                                             v = (-1)^(row$stimulus=="lower")*V[row$condition],
                                                             t0 = t0, z = z, sz = sz, st0=st0,
                                                             sv = SV[row$condition],
                                                             simult_conf = simult_conf,
                                                             z_absolute = FALSE, precision = precision))))


  p <- integrate(integrand, lower=paramDf$t0, upper=maxrt, subdivisions = subdivisions,
                 stop.on.error = stop.on.error)
  if (.progress) pb$tick()
  return(data.frame(p = p$value, info = p$message, err = p$abs.error))
}

