fill_thresholds <- function(res, used_cats, actual_nRatings, min_conf) {
  ### This function fills up the missing confidence thresholds with the best/easiest choices
  ### for the not identifiable or not reasonably fittable thresholds because some confidence
  ### categories are not used by the subject
  ### min_conf : theoretical minimum of confidence measure (depending on model: 0 or -Inf (then -1e+24 should be used))
  ### ToDo:   For sym_thetas==FALSE, use different nRatings for lower and upper responses in fitting
  symmetric_confidence_thresholds <- length(grep(pattern = "thetaUpper", names(res), value = T))<1
  if (symmetric_confidence_thresholds) {
    thetas <- rep(NA,(actual_nRatings-1))
    names(thetas) <- paste("theta", 1:(actual_nRatings-1), sep="")
    thetas[paste("theta", (used_cats[used_cats<max(used_cats)]), sep="")] <- c(t(res[,grep(pattern = "theta[0-9]", names(res), value = T)]))
    if (min(used_cats) >1) {
      thetas[paste("theta", 1:(min(used_cats)-1), sep="")] <- min(min_conf, min(thetas, na.rm = TRUE))
    }
    if (max(used_cats)<actual_nRatings) {
      thetas[paste("theta", (max(used_cats):(actual_nRatings-1)), sep="")] <- max(1e+24, max(thetas, na.rm=TRUE))
    }
    while (any(is.na(thetas))) {
      thetas[which(is.na(thetas))] <- thetas[(which(is.na(thetas))-1)]
    }
    res[,names(thetas)]<- thetas
  } else {
    thetasUpper <- rep(NA,(actual_nRatings-1))
    thetasLower <- rep(NA,(actual_nRatings-1))
    names(thetasUpper) <- paste("thetaUpper", 1:(actual_nRatings-1), sep="")
    names(thetasLower) <- paste("thetaLower", 1:(actual_nRatings-1), sep="")
    thetasUpper[paste("thetaUpper", (used_cats[used_cats<max(used_cats)]), sep="")] <- c(t(res[,grep(pattern = "thetaUpper", names(res), value = T)]))
    thetasLower[paste("thetaLower", (used_cats[used_cats<max(used_cats)]), sep="")] <- c(t(res[,grep(pattern = "thetaLower", names(res), value = T)]))
    if (min(used_cats) >1) {
      thetasUpper[paste("thetaUpper", 1:(min(used_cats)-1), sep="")] <- min(min_conf, min(thetasUpper, na.rm=TRUE))
      thetasLower[paste("thetaLower", 1:(min(used_cats)-1), sep="")] <- min(min_conf, min(thetasLower, na.rm=TRUE))
    }
    if (max(used_cats)<actual_nRatings) {
      thetasUpper[paste("thetaUpper", (max(used_cats):(actual_nRatings-1)), sep="")] <- max(1e+24, max(thetasUpper, na.rm = TRUE))
      thetasLower[paste("thetaLower", (max(used_cats):(actual_nRatings-1)), sep="")] <- max(1e+24, max(thetasLower, na.rm = TRUE))
    }
    while (any(is.na(thetasUpper))) {
      thetasUpper[which(is.na(thetasUpper))] <- thetasUpper[(which(is.na(thetasUpper))-1)]
    }
    while (any(is.na(thetasLower))) {
      thetasLower[which(is.na(thetasLower))] <- thetasLower[(which(is.na(thetasLower))-1)]
    }

    res[,names(thetasUpper)]<- thetasUpper
    res[,names(thetasLower)]<- thetasLower
  }
  res
}






### NICHT FERTIG!
fill_thresholdsDDM <- function(res, used_cats, actual_nRatings) {
  ### This function fills up the missing confidence thresholds with the best/easiest choices
  ### for the not identifiable or not reasonably fittable thresholds because some confidence
  ### categories are not used by the subject
  ### ToDo:   For sym_thetas==FALSE, use different nRatings for lower and upper responses in fitting
  symmetric_confidence_thresholds <- length(grep(pattern = "thetaUpper", names(res), value = T))<1
  if (symmetric_confidence_thresholds) {
    thetas <- rep(NA,(actual_nRatings-1))
    names(thetas) <- paste("theta", (actual_nRatings-1):1, sep="")
    thetas[paste("theta", sort((used_cats[used_cats<max(used_cats)]),decreasing=TRUE), sep="")] <- c(t(res[,grep(pattern = "theta[0-9]", names(res), value = T)]))
    if (min(used_cats) >1) {
      thetas[paste("theta", 1:(min(used_cats)-1), sep="")] <- max(1e+24, max(thetas, na.rm = TRUE))
    }
    if (max(used_cats)<actual_nRatings) {
      thetas[paste("theta", (max(used_cats):(actual_nRatings-1)), sep="")] <- min(0, min(thetas, na.rm=TRUE))
    }
    thetas[which(is.na(thetas))] <- thetas[(which(is.na(thetas))+1)]
    res[,names(thetas)]<- thetas
  } else {
    thetasUpper <- rep(NA,(actual_nRatings-1))
    thetasLower <- rep(NA,(actual_nRatings-1))
    names(thetasUpper) <- paste("thetaUpper", (actual_nRatings-1):1, sep="")
    names(thetasLower) <- paste("thetaLower", (actual_nRatings-1):1, sep="")
    thetasUpper[paste("thetaUpper", sort(used_cats[used_cats<max(used_cats)], decreasing = TRUE), sep="")] <- c(t(res[,grep(pattern = "thetaUpper", names(res), value = T)]))
    thetasLower[paste("thetaLower", sort(used_cats[used_cats<max(used_cats)], decreasing = TRUE), sep="")] <- c(t(res[,grep(pattern = "thetaLower", names(res), value = T)]))
    if (min(used_cats) >1) {
      thetasUpper[paste("thetaUpper", 1:(min(used_cats)-1), sep="")] <- max(1e+24, max(thetasUpper, na.rm=TRUE))
      thetasLower[paste("thetaLower", 1:(min(used_cats)-1), sep="")] <- max(1e+24, max(thetasLower, na.rm=TRUE))
    }
    if (max(used_cats)<actual_nRatings) {
      thetasUpper[paste("thetaUpper", (max(used_cats):(actual_nRatings-1)), sep="")] <- min(0, min(thetasUpper, na.rm = TRUE))
      thetasLower[paste("thetaLower", (max(used_cats):(actual_nRatings-1)), sep="")] <- min(0, min(thetasLower, na.rm = TRUE))
    }
    while (any(is.na(thetasUpper))) {
      thetasUpper[which(is.na(thetasUpper))] <- thetasUpper[(which(is.na(thetasUpper))+1)]
    }
    while (any(is.na(thetasLower))) {
      thetasLower[which(is.na(thetasLower))] <- thetasLower[(which(is.na(thetasLower))+1)]
    }

    res[,names(thetasUpper)]<- thetasUpper
    res[,names(thetasLower)]<- thetasLower
  }
  res
}

