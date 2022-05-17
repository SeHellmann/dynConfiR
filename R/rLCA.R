#' Simulation of confidence ratings and RTs in leaky competing accumulator model
#'
#' Simulates the decision responses, reaction times and state of the loosing accumulator
#' together with a discrete confidence judgment  in the leaky competing accumulator model.
#' Optionally, there is a post-decisional accumulation period, where the process continues.
#' Parameters required (in \code{paramDf}):
#'    mu1 and mu2 (mean momentary evidence for alternatives)
#'    pi (factor for input dependent noise of infinitesimal variance of processes)
#'    sig (input independent component of infinitesimal variance of processes)
#'    th (decision threshold)
#'    k (leakage)
#'    beta (inhibition)
#'    SPV (variation in starting points)
#'    tau (fixed post decisional accumulation period)
#'    wx (weight on balance of evidence in confidence measure)
#'    wrt (weight on RT in confidence measure)
#'    wint (weight on interaction of evidence and RT in confidence measure)
#'    t0 (minimal non-decision time)
#'    st0 (range of uniform distribution of non-decision time)
#' Also computes the Gamma rank correlation between the confidence
#' ratings and condition (task difficulty), reaction times and accuracy in the simulated output.
#'
#' @param paramDf a list or dataframe with one row. Column names should match the names of
#' model parameters (mu1, mu2, pi, sig, th, k, beta, SPV, and tau, wx, wrt, and wint and t0 and st0). For different
#' experimental conditions, simply numerate the respective parameters (e.g. for varying
#' drift rates, use mu11, mu12, mu13, mu21, mu22, mu23). All non-numerated parameters are
#' assumed to be constant across conditions. Additionally, the confidence thresholds may be given by names with
#' thetaUpper1, thetaUpper2,..., thetaLower1,... or, for symmetric thresholds only by theta1, theta2,....
#' Note that confidence thresholds are not allowed to vary with experimental condition.
#' @param n integer. The number of samples (per condition (and stimulus direction)) generated.
#' Total number of samples is \code{n*nConditions*length(stimulus)}.
#' @param time_scaled logical. Whether a time_scaled transformation for the confidence measure should
#' be used.
#' @param gamma logical. If TRUE, the gamma correlation between confidence ratings, rt and accuracy is
#' computed.
#' @param agg_simus logical. Simulation is done on a trial basis with rts outcome. If TRUE,
#' the simulations will be aggregated over RTs to return only the distribution of response and
#' confidence ratings. Default: FALSE.
#' @param stimulus numeric vector. Either 1, 2 or c(1, 2) (default).
#' Together with condition represents the experimental situation. In a 2AFC task the presented
#' stimulus belongs to one of two categories. In the default setting trials with
#' both categories presented are simulated but one can choose to simulate only trials with the
#' stimulus coming from one category (each associated with positive drift in one of two accumulators).
#' More precisely, if stimulus is 2, trials with interchanged mu1 and mu2 are simulated representing a trial
#' with contrary stimulus identity as represented in the paramDf.
#' @param delta numerical. Size of steps for the discretized simulation (see details).
#' @param maxrt numerical. Maximum reaction time to be simulated (see details). Default: 15.
#' @param seed numerical. Seeding for non-random data generation. (Also possible outside of the function.)
#'
#' @return Depending on gamma and agg_simus. If gamma is TRUE, returns a list with elements:
#' "simus" (the simulated data frame) and "gamma", which is again a list with elements
#' "condition", "rt" and "correct", each a tibble with two columns (see details for more
#' information). If gamma is FALSE, returns a data frame with columns: condition, stimulus,
#' response, correct, rt, conf (the continuous confidence measure) and rating (the discrete
#' confidence rating) or (if agg_simus=TRUE): condition, stimulus, response, correct, rating
#' and p (for the probability of a response and rating, given the condition and stimulus).
#'
#'
#' @details The simulation is done by simulating discretized steps until one process reaches
#' the boundary with an update rule:
#' \deqn{\delta X_i(t) = \max (0, X_i(t) + \delta_t ((k-1)X_i(t)-\beta X_{j=i} (t) + \mu_i + \varepsilon_i (t)),}
#' with \eqn{\varepsilon_i(t) \sim N(0, (\pi \mu_i)^2 + \sigma^2 )}. If no boundary is met within the maximum time, response is
#' set to 0. After the decision, the accumulation continues for a time period (tau), until
#' the final state is used for the computation of confidence.
#'
#' The Gamma coefficients are computed separately for
#' correct/incorrect responses for the correlation of confidence ratings with condition and rt
#' and separately for conditions for the correlation of accuracy and confidence. The resulting
#' tibbles in the output thus have two columns. One for the grouping variable and one for the
#' Gamma coefficient.
#'
#' @note Different parameters for different conditions are only allowed for drift rate, \code{v},
#' and variability, \code{s}. All other parameters are used for all conditions.
#'
#'
#'
#' @author Sebastian Hellmann.
#'
#' @name rLCA
#' @import dplyr
#' @importFrom magrittr %>%
#' @importFrom Hmisc rcorr.cens
#' @importFrom rlang .data
#' @importFrom stats runif
# @importFrom pracma integral
#' @aliases simulateLCA
#'

## When given vectorised parameters, n is the number of replicates for each parameter set
#' @rdname rLCA
#' @export
rLCA <- function (paramDf, n=1e+4,   time_scaled=FALSE,
                        gamma = FALSE, agg_simus=FALSE,
                        stimulus = c(1,2), delta=0.01, maxrt=15, seed=NULL)
{
  # library(tidyverse)
  # n <- 100
  # time_scaled = FALSE
  # gamma = FALSE
  # agg_simus = FALSE
  # stimulus = c(1,2)
  # delta = 0.01
  # maxrt = 15
  # seed = NULL
  #
  # paramDf <- as.data.frame(list("mu1"=1, "mu2"=0, "pi"=1, "sig"=1,
  #                 "th"=2, "k"=0.9, "beta"=0.2,
  #                 "SPV"=0.1, "tau"=1,
  #                 "theta1" = 0.2, "theta2" = 0.4, "theta3" = 1,
  #                 "wx" = 0.7, "wrt" = 0.2, "wint"=0.1,
  #                   "t0" = 0.1, "st0"=0.05))
  #
  #
  # paramDf <- as.data.frame(list("mu11"=1, "mu21"=0,
  #                 "mu12"=2, "mu22"=1,
  #                 "pi"=1, "sig"=1,
  #                 "th"=2, "k"=0.9, "beta"=0.2,
  #                 "SPV"=0.1, "tau"=1,
  #                 "theta1" = 0.2, "theta2" = 0.4, "theta3" = 1,
  #                 "wx" = 0.7, "wrt" = 0.2, "wint"=0.1,
  #                   "t0" = 0.1, "st0"=0.05))

  ## Create dummy parameters, s.t. R CMD check does not complain that they are not defined
  mu1 <- 0
  mu2 <- 0
  sig <- 0
  th <- 0
  k <- 0
  SPV <- 0
  tau <- 0
  wx <- 0
  wint <- 0
  wrt <- 0

  if (!is.null(seed)) {
    set.seed(seed)
  }

  if (!(all(stimulus %in% c(1,2)))) {
    stop(paste("Not accepted value for stimulus: ", paste(stimulus, collapse=", "),". Must be either 1, 2 or c(1,2).", sep=""))
  }
  if (("wx1" %in% names(paramDf)) || ("wrt1" %in% names(paramDf)) || ("wint1" %in% names(paramDf))) {
    stop("Weights in confidence measure should not change with conditions (at least its not implemented, yet)!")
  }
  if (!time_scaled) {
    paramDf$wint <- 0
    paramDf$wrt <- 0
    paramDf$wx <- 1
  }

  ## recover parameters from paramDf
  bin_conf = FALSE
  if (length(grep(pattern = "theta", names(paramDf)))>=1) {
    bin_conf = TRUE
    symmetric_confidence_thresholds <- length(grep(pattern = "thetaUpper", names(paramDf), value = T))<1
    if (symmetric_confidence_thresholds) {
      nRatings <- length(grep(pattern = "^theta[0-9]", names(paramDf)))+1
      thetas_1 <- c(-Inf, t(paramDf[paste("theta",1:(nRatings-1), sep = "")]), Inf)
      thetas_2 <- c(-Inf, t(paramDf[paste("theta",1:(nRatings-1), sep = "")]), Inf)
    } else {
      nRatings <- length(grep(pattern = "^thetaUpper[0-9]", names(paramDf)))+1
      thetas_1 <- c(-Inf, t(paramDf[paste("thetaUpper",1:(nRatings-1), sep = "")]), Inf)
      thetas_2 <- c(-Inf, t(paramDf[paste("thetaLower",1:(nRatings-1), sep="")]), Inf)
    }
    paramDf <- paramDf[names(paramDf)[!grepl("theta", names(paramDf))]]
  }
  parnames <- c("mu1", "mu2", "pi", "sig", "th","k", "beta", "SPV", "tau", "wx", "wrt", "wint", "t0", "st0")
  cond_pars <- NULL
  for (l in 1:length(parnames)) {
    if (length(grep(pattern=paste0(parnames[l], "[0-9]"), names(paramDf)))>1) {
      cond_pars <- c(cond_pars, parnames[l])
    }
  }
  if (length(cond_pars) > 0 ) {
    nConds <- length(grep(pattern = paste0(cond_pars[1],"[0-9]"), names(paramDf), value = T))
    df <- expand.grid(condition = 1:nConds, stimulus=stimulus)
    for (l in 1:length(cond_pars)) {
      #df[[cond_pars[l]]] <- NA
      for (i in 1:nrow(df)) {
        df[[cond_pars[l]]][i] <- paramDf[[paste0(cond_pars[l], df[i, "condition"])]]
      }
      paramDf <- paramDf[!grepl(pattern=cond_pars[l], names(paramDf))]
      #df[[cond_pars[l]]] <- paramDf[grep(pattern=cond_pars[l], names(paramDf), value = TRUE)][df$condition]
    }
  } else {
    nConds <- 1
  }
  if (!("mu1" %in% names(df))) {
    df$mu1 <- paramDf["mu1"]
    df$mu2 <- paramDf["mu2"]
  }
  if (length(stimulus) > 1) {
    temp <- df[df$stimulus==2, "mu1"]
    df[df$stimulus==2, "mu1"] <- df[df$stimulus==2, "mu2"]
    df[df$stimulus==2, "mu2"] <- temp
  }
  df$stimulus <- if_else(df$mu1==df$mu2, 0,
                         if_else(df$mu1>df$mu2, 1, 2))
  df <- df[!duplicated(df),]

  paramDf[c("wx", "wrt", "wint")] <-   paramDf[c("wx", "wrt", "wint")]/sum(as.numeric(  paramDf[c("wx", "wrt", "wint")]))

  #list2env(paramDf, envir = environment())
  parnames <- names(paramDf)
  for (i in 1:length(parnames)) {
    assign(parnames[i], paramDf[[parnames[i]]], envir = environment())
  }

  help_fct <- function(mu1, mu2, sig, pi, th, k, beta, SPV, tau) {
    res <- as.data.frame(r_LCA(n,
                               c(mu1, mu2, sig, pi, th, k, beta, SPV, tau),
                              delta=delta, maxT=maxrt))
    names(res) <- c("rt", "response", "xl", "x1", "x2")
    res
  }
  ## Produce process outcomes and compute confidence measure
  simus <- df %>%
    group_by(.data$condition, .data$stimulus) %>%
    summarise(help_fct(mu1, mu2, sig, pi, th, k, beta, SPV, tau)) %>%
    mutate(BoE = (.data$x2-.data$x1)*(if_else(.data$response==2, 1, -1)),
           conf = if_else(rep(time_scaled, n),
                          wx*.data$BoE + wrt/sqrt(.data$rt) + wint*.data$BoE/sqrt(.data$rt),
                          wx*.data$BoE))
  if ("t0" %in% cond_pars) {
    t0 <- df[1:nConds,"t0"][simus$condition]
  }
  if ("st0" %in% cond_pars) {
    st0 <- df[1:nConds,"st0"][simus$condition]
  }
  simus <- simus %>% ungroup() %>%
    rename(dt = "rt") %>%
    mutate(rt = .data$dt + runif(nrow(simus), min = t0, max=t0+st0),
           correct = if_else(.data$stimulus==0, 3, as.numeric(.data$response==.data$stimulus)))

  ### Bin confidence measure for discrete ratings, if parameters given:
  if (!bin_conf) {
    if (gamma) {
      warning("gamma correlations only returned, if theta-parameters given and confidence returned")
    }
    if (agg_simus) {
      simus <- simus %>% group_by(.data$correct, .data$condition) %>%
        summarise(p = n()/(2*n)) %>%
        full_join(y=expand.grid(condition=1:nConds,
                                correct=c(0,1))) %>%
        mutate(p = ifelse(is.na(.data$p), 0, .data$p))
    } else {
      simus <- simus[c("condition", "stimulus", "response", "correct", "rt","BoE", "conf", "dt")]
    }
    return(simus)
  }


  simus$rating <- 1
  simus$rating[simus$response==1] <- as.numeric(as.factor(cut(simus$conf[simus$response==1], breaks=thetas_1)))
  simus$rating[simus$response==2] <- as.numeric(as.factor(cut(simus$conf[simus$response==2], breaks=thetas_2)))




  if (gamma==TRUE) {
    gamma_condition <- simus %>% group_by(.data$correct) %>%
      summarise(data.frame(t(rcorr.cens(.data$rating,S=.data$condition, outx=TRUE)))) %>%
      select(.data$correct, Gamma = .data$Dxy)
    gamma_rt <- simus %>% group_by(.data$correct) %>%
      summarise(data.frame(t(rcorr.cens(.data$rating,S=.data$rt, outx=TRUE))))%>%
      select(.data$correct, Gamma = .data$Dxy)
    gamma_correct <- simus %>% group_by(.data$condition) %>%
      summarise(data.frame(t(rcorr.cens(.data$rating,S=.data$correct, outx=TRUE))))%>%
      select(.data$condition, Gamma = .data$Dxy)
    gamma_rt_bycondition <- simus %>% group_by(.data$condition) %>%
      summarise(data.frame(t(rcorr.cens(.data$rating,S=.data$rt, outx=TRUE))))%>%
      select(.data$condition, Gamma = .data$Dxy)
    gamma_rt_byconditionbycorrect <- simus %>% group_by(.data$condition, .data$correct) %>%
      summarise(data.frame(t(rcorr.cens(.data$rating,S=.data$rt, outx=TRUE))))%>%
      select(.data$condition, Gamma = .data$Dxy)
  }
  if (agg_simus) {
    simus <- simus %>% group_by(.data$rating, .data$correct, .data$condition) %>%
      summarise(p = n()/(2*n)) %>%
      full_join(y=expand.grid(rating=1:nRatings, condition=1:nConds,
                              correct=c(0,1))) %>%
      mutate(p = ifelse(is.na(.data$p), 0, .data$p))
  } else {
    simus <- simus[c("condition", "stimulus", "response", "rating", "correct", "rt","BoE", "conf", "dt")]
  }
  if (gamma) {
    return(list("simus"=simus,
                "gamma" = list("condition" = gamma_condition,
                               "rt" = gamma_rt,
                               "correct" = gamma_correct,
                               "rt_bycondition" = gamma_rt_bycondition,
                               "rt_byconditionbycorrect" = gamma_rt_byconditionbycorrect)))
  } else {
    return(simus)
  }
}


