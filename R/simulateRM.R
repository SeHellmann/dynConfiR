#' Simulation of confidence ratings and RTs in race confidence models
#'
#' Simulates the decision responses, reaction times and state of the loosing accumulator
#' together with a discrete confidence judgment  in the independent and partially anti-correlated
#' race model (IRM and PCRM) (Hellmann et al., 2023), given specific parameter constellations.
#' See \link{RaceModels} for more information about
#' parameters. Also computes the Gamma rank correlation between the confidence
#' ratings and condition (task difficulty), reaction times and accuracy in the
#' simulated output. Basically, this function is a wrapper for \code{\link{rIRM}}
#' and \code{\link{rPCRM}} for application in confidence experiments with
#' manipulation of specific parameters.
#' `rRM_Kiani` simulates a different version of race models, presented in
#' Kiani et al. (2014), but without a confidence measure.
#'
#' @param paramDf a list or data frame with one row. Column names should match the names of
#' \link{RaceModels} parameter names (only `mu1` and `mu2` are not used in this context but
#' replaced by the parameter `v`). For different stimulus quality/mean
#' drift rates, names should be `v1`, `v2`, `v3`,....
#' Different `s` parameters are possible with `s1`, `s2`, `s3`,... with equally many steps as for drift rates. Additionally, the confidence
#' thresholds should be given by names with `thetaUpper1`, `thetaUpper2`,..., `thetaLower1`,... or,
#' for symmetric thresholds only by `theta1`, `theta2`,....
#' @param n integer. The number of samples (per condition and stimulus direction) generated.
#' Total number of samples is \code{n*nConditions*length(stimulus)}.
#' @param model character scalar. One of "IRM" or "PCRM". ("IRMt" and "PCRMt" will also be accepted. In that case,
#' time_scaled is set to TRUE.)
#' @param time_scaled logical. Whether a time_scaled transformation for the confidence measure should
#' be used.
#' @param gamma logical. If TRUE, the gamma correlation between confidence ratings, rt and accuracy is
#' computed.
#' @param agg_simus logical. Simulation is done on a trial basis with RTs outcome. If TRUE,
#' the simulations will be aggregated over RTs to return only the distribution of response and
#' confidence ratings. Default: FALSE.
#' @param stimulus numeric vector. Either 1, 2 or c(1, 2) (default).
#' Together with condition represents the experimental situation. In a binary decision task the presented
#' stimulus belongs to one of two categories. In the default setting trials with
#' both categories presented are simulated but one can choose to simulate only trials with the
#' stimulus coming from one category (each associated with positive drift in one of two accumulators).
#' @param delta numerical. Size of steps for the discretized simulation (see details).
#' @param maxrt numerical. Maximum reaction time to be simulated (see details). Default: 15.
#' @param seed numerical. Seeding for non-random data generation. (Also possible outside of the function.)
#'
#' @return Depending on `gamma` and `agg_simus`.
#'
#' If `gamma` is `FALSE`, returns a `data.frame` with columns: `condition`,
#' `stimulus`, `response`, `correct`, `rt`, `conf` (the continuous confidence
#' measure) and `rating` (the discrete confidence rating) or
#' (if `agg_simus=TRUE`): `condition`, `stimulus`,`response`, `correct`,
#' `rating` and `p` (for the probability of a response and rating, given
#' the condition and stimulus).
#'
#' If `gamma` is `TRUE`, returns a `list` with elements:
#' `simus` (the simulated data frame) and `gamma`, which is again a `list` with elements
#' `condition`, `rt` and `correct`, each a `tibble` with two columns (see details for more
#' information).
#'
#'
#' @details The simulation is done by simulating normal variables in discretized steps until
#' one process reaches the boundary. If no boundary is met within the maximum time, response is
#' set to 0. The output of the fitting function \code{\link{fitRTConf}} with the respective model
#' fits the argument `paramDf` for simulation. The Gamma coefficients are computed separately for
#' correct/incorrect responses for the correlation of confidence ratings with condition and rt
#' and separately for conditions for the correlation of accuracy and confidence. The resulting
#' data frames in the output thus have two columns. One for the grouping variable and one for the
#' Gamma coefficient.
#'
#' @note Different parameters for different conditions are only allowed for drift rate, \code{v},
#' and process variability, \code{s}. All other parameters are used for all conditions.
#'
#' @references  Hellmann, S., Zehetleitner, M., & Rausch, M. (2023). Simultaneous modeling of choice, confidence and response time in visual perception. \emph{Psychological Review} 2023 Mar 13. doi: 10.1037/rev0000411. Epub ahead of print. PMID: 36913292.
#'
#' Kiani, R., Corthell, L., & Shadlen, M.N. (2014) Choice certainty is informed
#' by both evidence and decision time.
#' Neuron, 84(6), 1329-1342. doi:10.1016/j.neuron.2014.12.015
#'
#' @author Sebastian Hellmann.
#'
#' @name simulateRM
#' @import dplyr
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @importFrom stats runif
#' @importFrom Rcpp evalCpp
# @importFrom pracma integral
#' @aliases simulateIRM simulatePCRM
#'
#' @examples
#' # Examples for "PCRM" model (equivalent applicable for "IRM" model)
#' # 1. Define some parameter set in a data.frame
#' paramDf <- data.frame(a=2,b=2, v1=0.5, v2=1, t0=0.1,st0=0,
#'                       wx=0.6, wint=0.2, wrt=0.2,
#'                       theta1=4)
#'
#' # 2. Simulate trials for both stimulus categories and all conditions (2)
#' simus <- simulateRM(paramDf, n=30,model="PCRM", time_scaled=TRUE)
#' head(simus)
#' # equivalent:
#' simus <- simulateRM(paramDf, model="PCRMt")
#' \donttest{
#'   library(ggplot2)
#'   simus <- simus[simus$response!=0,]
#'   simus$rating <- factor(simus$rating, labels=c("unsure", "sure"))
#'   ggplot(simus, aes(x=rt, group=interaction(correct, rating),
#'                     color=as.factor(correct), linetype=rating))+
#'     geom_density(size=1.2)+
#'     facet_grid(rows=vars(condition), labeller = "label_both")
#' }
#'
#' # automatically aggregate simulation distribution
#' # to get only accuracy x confidence rating distribution for
#' # all conditions
#' agg_simus <- simulateRM(paramDf, n = 20, model="PCRMt", agg_simus = TRUE)
#' head(agg_simus)
#' \donttest{
#'   agg_simus$rating <- factor(agg_simus$rating, labels=c("unsure", "sure"))
#'   library(ggplot2)
#'   ggplot(agg_simus, aes(x=rating, group=correct, fill=as.factor(correct), y=p))+
#'     geom_bar(stat="identity", position="dodge")+
#'     facet_grid(cols=vars(condition), labeller = "label_both")
#' }
#'
#' \donttest{
#'   # Compute Gamma correlation coefficients between
#'   # confidence and other behavioral measures
#'   # output will be a list
#'   simu_list <- simulateRM(paramDf, model="IRMt", gamma=TRUE)
#'   simu_list
#' }
#'


## When given vectorised parameters, n is the number of replicates for each parameter set
#' @rdname simulateRM
#' @export
simulateRM <- function (paramDf, n=1e+4,  model = "IRM", time_scaled=FALSE,
                        gamma = FALSE, agg_simus=FALSE,
                        stimulus = c(1,2), delta=0.01, maxrt=15, seed=NULL)
{
  if (gamma && !requireNamespace("Hmisc", quietly = TRUE)) {
    warning("Package 'Hmisc' is not installed, but required to computed Gamma correlations.
    Please install 'Hmisc', if computation of Gamma is required.
    Otherwise the function will continue without computation of Gamma.
    (Computation of Gamma is still possible with the function output, if agg_simus = FALSE.)",
    immediate.=TRUE)
    gamma <- FALSE
  }
  if (!is.null(seed)) {
    set.seed(seed)
  }
  #### Check model argument
  if (model=="IRMt") {
    model = "IRM"
    time_scaled=TRUE
  }
  if (model=="PCRMt") {
    model = "PCRM"
    time_scaled=TRUE
  }
  if (!model %in% c("IRM", "PCRM")) stop("model must be 'IRM', 'PCRM', 'IRMt' or 'PCRMt'")

  if (!(all(stimulus %in% c(1,2)))) {
    stop(paste("Not accepted value for stimulus: ", paste(stimulus, collapse=", "),". Must be either 1, 2 or c(1,2).", sep=""))
  }


  ## recover parameters from paramDf

  nConds <- length(grep(pattern = "^v[0-9]", names(paramDf), value = T))
  if (nConds > 0 ) {
    V <- c(t(paramDf[,paste("v",1:(nConds), sep = "")]))
  } else {
    V <- paramDf$v
    nConds <- 1
  }
  vary_s <-   length(grep(pattern = "^s[0-9]", names(paramDf), value = T))>1
  if (vary_s){
    S <- c(t((paramDf[,paste("s",1:(nConds), sep = "")])))
  } else {
    if ("s" %in% names(paramDf)) {
      S <- rep(paramDf$s, nConds)
    } else {
      S <- rep(1, nConds)
    }
  }
  df <- expand.grid(condition = 1:nConds, stimulus=stimulus)
  ## Produce process outcomes and compute confidence measure
  simus <- data.frame()
  for ( i in 1:nrow(df)) {
    mu1 <- (-1)^(1+as.numeric(df[i,]$stimulus==1)) * V[df[i,]$condition]
    temp <- as.data.frame(r_RM(n,c(mu1, -mu1,
                   -paramDf$a, -paramDf$b,
                   S[df[i,]$condition],S[df[i,]$condition],
                   0,0,0,0),
                   rho = ifelse(model=="IRM", 0, -.5),
                   delta=delta, maxT=maxrt))
    names(temp) <- c("rt", "response", "xl")
    simus <- rbind(simus,
                 cbind(condition=df[i, "condition"], stimulus=df[i, "stimulus"], temp))
  }
  simus$xj <- if_else(simus$response==1, paramDf$b, paramDf$a) + simus$xl
  if (time_scaled) {
    simus$conf <- -paramDf$wx*simus$xl + (paramDf$wrt / sqrt(simus$rt)) - paramDf$wint * (simus$xl / sqrt(simus$rt))
  } else {
    simus$conf <- -simus$xl
  }
  simus$rt <- simus$rt + runif(nrow(simus), min = paramDf$t0, max=paramDf$t0+paramDf$st0)
  simus$correct <- as.numeric(simus$response == simus$stimulus)




  ### Bin confidence measure for discrete ratings, if parameters given:
  if (length(grep(pattern = "theta", names(paramDf)))<1) {
    if (gamma) {
      warning("gamma correlations only returned, if theta-parameters given and confidence returned")
    }
    if (agg_simus) {
      simus <- simus %>% group_by(.data$correct, .data$condition) %>%
        summarise(p = n()/(2*n), .groups = "drop") %>%
        full_join(y=expand.grid(condition=1:nConds,
                                correct=c(0,1)), by=join_by("correct","condition")) %>%
        mutate(p = ifelse(is.na(.data$p), 0, .data$p))
    } else {
      simus <- simus[c("condition", "stimulus", "response", "correct", "rt","xj")]
    }
    return(simus)
  }

  symmetric_confidence_thresholds <- length(grep(pattern = "thetaUpper", names(paramDf), value = T))<1
  if (symmetric_confidence_thresholds) {
    nRatings <- length(grep(pattern = "^theta[0-9]", names(paramDf)))+1
    if (nRatings==1) nRatings <- 2
    thetas_upper <- c(-Inf, t(paramDf[,paste("theta",1:(nRatings-1), sep = "")]), Inf)
    thetas_lower <- c(-Inf, t(paramDf[,paste("theta",1:(nRatings-1), sep = "")]), Inf)
  } else {
    nRatings <- length(grep(pattern = "^thetaUpper[0-9]", names(paramDf)))+1
    if (nRatings==1) nRatings <- 2
    thetas_upper <- c(-Inf, t(paramDf[,paste("thetaUpper",1:(nRatings-1), sep = "")]), Inf)
    thetas_lower <- c(-Inf, t(paramDf[,paste("thetaLower",1:(nRatings-1), sep="")]), Inf)
  }
  levels_lower <- cumsum(as.numeric(table(thetas_lower)))
  levels_lower <- levels_lower[-length(levels_lower)]
  levels_upper <- cumsum(as.numeric(table(thetas_upper)))
  levels_upper <- levels_upper[-length(levels_upper)]
  thetas_lower <- unique(thetas_lower)
  thetas_upper <- unique(thetas_upper)

  simus$rating <- 1
  simus$rating[simus$response==1] <- as.numeric(as.character(cut(simus$conf[simus$response==1],
                                                                        breaks=thetas_upper, labels = levels_upper)))
  simus$rating[simus$response==2] <- as.numeric(as.character(cut(simus$conf[simus$response==2],
                                                                        breaks=thetas_lower, labels = levels_lower)))




  if (gamma==TRUE) {
    gamma_condition <- simus %>% group_by(.data$correct) %>%
      summarise(Gamma=c(t(Hmisc::rcorr.cens(.data$rating,S=.data$condition, outx=TRUE)))[2])
    gamma_rt <- simus %>% group_by(.data$correct) %>%
      summarise(Gamma=c(t(Hmisc::rcorr.cens(.data$rating,S=.data$rt, outx=TRUE)))[2])
    gamma_correct <- simus %>% group_by(.data$condition) %>%
      summarise(Gamma=c(t(Hmisc::rcorr.cens(.data$rating,S=.data$correct, outx=TRUE)))[2])
    gamma_rt_bycondition <- simus %>% group_by(.data$condition) %>%
      summarise(Gamma=c(t(Hmisc::rcorr.cens(.data$rating,S=.data$rt, outx=TRUE)))[2])
    gamma_rt_byconditionbycorrect <- simus %>% group_by(.data$condition, .data$correct) %>%
      summarise(Gamma=c(t(Hmisc::rcorr.cens(.data$rating,S=.data$rt, outx=TRUE)))[2])
  }
  if (agg_simus) {
    simus <- simus %>% group_by(.data$rating, .data$correct, .data$condition) %>%
      summarise(p = n()/(2*n)) %>%
      full_join(y=expand.grid(rating=1:nRatings, condition=1:nConds,
                              correct=c(0,1)), by=join_by("rating", "condition", "correct")) %>%
      mutate(p = ifelse(is.na(.data$p), 0, .data$p))
  } else {
    simus <- simus[c("condition", "stimulus", "response", "correct", "rt","xj",  "conf", "rating")]
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


## When given vectorised parameters, n is the number of replicates for each parameter set
#' @rdname simulateRM
#' @export
rRM_Kiani <- function (paramDf, n=1e+4, time_scaled=FALSE,
                 gamma = FALSE, agg_simus=FALSE,
                 stimulus = c(1,2), delta=0.01, maxrt=15, seed=NULL)
{
  if (gamma && !requireNamespace("Hmisc", quietly = TRUE)) {
    warning("Package 'Hmisc' is not installed, but required to computed Gamma correlations.
    Please install 'Hmisc', if computation of Gamma is required.
    Otherwise the function will continue without computation of Gamma.
    (Computation of Gamma is still possible with the function output, if agg_simus = FALSE.)",
            immediate.=TRUE)
    gamma <- FALSE
  }
  if (!is.null(seed)) {
    set.seed(seed)
  }
  if (!("rho" %in% names(paramDf) && "Bl" %in% names(paramDf))) {
    stop("Parameters rho and Bl must be supplied for Kiani's model. At least one is missing.")
  }
  if (!(all(stimulus %in% c(1,2)))) {
    stop(paste("Not accepted value for stimulus: ", paste(stimulus, collapse=", "),". Must be either 1, 2 or c(1,2).", sep=""))
  }


  ## recover parameters from paramDf

  nConds <- length(grep(pattern = "^v[0-9]", names(paramDf), value = T))
  if (nConds > 0 ) {
    V <- c(t(paramDf[,paste("v",1:(nConds), sep = "")]))
  } else {
    V <- paramDf$v
    nConds <- 1
  }
  vary_s <-   length(grep(pattern = "^s[0-9]", names(paramDf), value = T))>1
  if (vary_s){
    S <- c(t((paramDf[,paste("s",1:(nConds), sep = "")])))
  } else {
    S <- rep(paramDf$s, nConds)
  }

  df <- expand.grid(condition = 1:nConds, stimulus=stimulus)
  ## Produce process outcomes and compute confidence measure
  simus <- data.frame()
  for ( i in 1:nrow(df)) {
    mu1 <- (-1)^(1+as.numeric(df[i,]$stimulus==1)) * V[df[i,]$condition]
    temp <- as.data.frame(r_RM_Kiani(n,c(mu1, -mu1, paramDf$a, paramDf$b, paramDf$s, paramDf$t0, paramDf$st0),
                                     paramDf$rho, paramDf$Bl,
                                     delta=delta, maxT=maxrt))
    names(temp) <- c("rt", "response", "xl")
    simus <- rbind(simus,
                   cbind(condition=df[i, "condition"], stimulus=df[i, "stimulus"], temp))
  }
  simus$xj <-simus$xl
  if (time_scaled) {
    simus$conf <- -paramDf$wx*simus$xl + (paramDf$wrt / sqrt(simus$rt)) - paramDf$wint * (simus$xl / sqrt(simus$rt))
  } else {
    simus$conf <- -simus$xl
  }

  simus$rt <- simus$rt + runif(nrow(simus), min = paramDf$t0, max=paramDf$t0+paramDf$st0)
  simus$correct <- as.numeric(simus$response == simus$stimulus)




  ### Bin confidence measure for discrete ratings, if parameters given:
  if (length(grep(pattern = "theta", names(paramDf)))<1) {
    if (gamma) {
      warning("gamma correlations only returned, if theta-parameters given and confidence returned")
    }
    if (agg_simus) {
      simus <- simus %>% group_by(.data$correct, .data$condition) %>%
        summarise(p = n()/(2*n)) %>%
        full_join(y=expand.grid(condition=1:nConds,
                                correct=c(0,1)), by=join_by("condition", "correct")) %>%
        mutate(p = ifelse(is.na(.data$p), 0, .data$p))
    } else {
      simus <- simus[c("condition", "stimulus", "response", "correct", "rt","xj")]
    }
    return(simus)
  }

  symmetric_confidence_thresholds <- length(grep(pattern = "thetaUpper", names(paramDf), value = T))<1
  if (symmetric_confidence_thresholds) {
    nRatings <- length(grep(pattern = "^theta[0-9]", names(paramDf)))+1
    thetas_upper <- c(-Inf, t(paramDf[,paste("theta",1:(nRatings-1), sep = "")]), Inf)
    thetas_lower <- c(-Inf, t(paramDf[,paste("theta",1:(nRatings-1), sep = "")]), Inf)
  } else {
    nRatings <- length(grep(pattern = "^thetaUpper[0-9]", names(paramDf)))+1
    thetas_upper <- c(-Inf, t(paramDf[,paste("thetaUpper",1:(nRatings-1), sep = "")]), Inf)
    thetas_lower <- c(-Inf, t(paramDf[,paste("thetaLower",1:(nRatings-1), sep="")]), Inf)
  }
  levels_lower <- cumsum(as.numeric(table(thetas_lower)))
  levels_lower <- levels_lower[-length(levels_lower)]
  levels_upper <- cumsum(as.numeric(table(thetas_upper)))
  levels_upper <- levels_upper[-length(levels_upper)]
  thetas_lower <- unique(thetas_lower)
  thetas_upper <- unique(thetas_upper)

  simus$rating <- 1
  simus$rating[simus$response==1] <- as.numeric(as.character(cut(simus$conf[simus$response==1],
                                                                        breaks=thetas_upper, labels = levels_upper)))
  simus$rating[simus$response==2] <- as.numeric(as.character(cut(simus$conf[simus$response==2],
                                                                        breaks=thetas_lower, labels = levels_lower)))



  if (gamma==TRUE) {
    gamma_condition <- simus %>% group_by(.data$correct) %>%
      summarise(Gamma=c(t(Hmisc::rcorr.cens(.data$rating,S=.data$condition, outx=TRUE)))[2])
    gamma_rt <- simus %>% group_by(.data$correct) %>%
      summarise(Gamma=c(t(Hmisc::rcorr.cens(.data$rating,S=.data$rt, outx=TRUE)))[2])
    gamma_correct <- simus %>% group_by(.data$condition) %>%
      summarise(Gamma=c(t(Hmisc::rcorr.cens(.data$rating,S=.data$correct, outx=TRUE)))[2])
    gamma_rt_bycondition <- simus %>% group_by(.data$condition) %>%
      summarise(Gamma=c(t(Hmisc::rcorr.cens(.data$rating,S=.data$rt, outx=TRUE)))[2])
    gamma_rt_byconditionbycorrect <- simus %>% group_by(.data$condition, .data$correct) %>%
      summarise(Gamma=c(t(Hmisc::rcorr.cens(.data$rating,S=.data$rt, outx=TRUE)))[2])
  }
  if (agg_simus) {
    simus <- simus %>% group_by(.data$rating, .data$correct, .data$condition) %>%
      summarise(p = n()/(2*n)) %>%
      full_join(y=expand.grid(rating=1:nRatings, condition=1:nConds,
                              correct=c(0,1)), by=join_by("rating", "condition", "correct")) %>%
      mutate(p = ifelse(is.na(.data$p), 0, .data$p))
  } else {
    simus <- simus[c("condition", "stimulus", "response", "correct", "rt","xj",  "conf", "rating")]
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





