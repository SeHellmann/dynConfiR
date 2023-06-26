#' Simulation of confidence ratings and RTs in dynWEV and 2DSD confidence models
#'
#' Simulates the decision responses and reaction times together with a
#' discrete confidence judgment in the dynaViTE model, the 2DSD model (Pleskac & Busemeyer, 2010)
#' and the dynWEV model (Hellmann et al., 2023), given specific parameter constellations.
#' See \code{\link{dWEV}} and \code{\link{d2DSD}} for more information about parameters.
#' Also computes the Gamma rank correlation between the confidence ratings and condition
#' (task difficulty), reaction times and accuracy in the simulated output.
#' Basically, this function is a wrapper for \code{\link{rWEV}} and \code{\link{r2DSD}}
#' for application in confidence experiments with manipulation of specific parameters.
#'
#' @param paramDf a list or dataframe with one row. Column names should match the names
#' of \link{dynaViTE} and \link{2DSD} model specific parameter names. For different stimulus quality/mean
#' drift rates, names should be `v1`, `v2`, `v3`,....
#' Different `sv` and/or `s` parameters are possible with `sv1`, `sv2`, `sv3`... (`s1`, `s2`, `s3`,...
#' respectively) with equally many steps as for drift rates. Additionally, the confidence
#' thresholds should be given by names with `thetaUpper1`, `thetaUpper2`,..., `thetaLower1`,... or,
#' for symmetric thresholds only by `theta1`, `theta2`,....
#' @param n integer. The number of samples (per condition and stimulus direction) generated.
#' Total number of samples is \code{n*nConditions*length(stimulus)}.
#' @param model character scalar. One of "dynaViTE", "dynWEV", or "2DSD".
#' @param delta numeric. Discretization steps for simulations with the stochastic process.
#' @param maxrt numeric. Maximum reaction time returned.
#' If the simulation of the stochastic process exceeds a rt of `maxrt`,
#' the response will be set to 0 and `maxrt` will be returned as rt.
#' @param  simult_conf logical. `TRUE`, if in the experiment confidence was reported simultaneously
#' with the decision, as then decision and confidence judgment are assumed to have happened
#' subsequent before response and tau is added to the simulated decision time. If `FALSE`
#' returned response time will only be decision time plus non-judgment time component.
#' @param gamma logical. If TRUE, the gamma correlation between confidence ratings, rt
#' and accuracy is computed.
#'
#' @param agg_simus logical. Simulation is done on a trial basis with RTs outcome.
#' If TRUE, the simulations will be aggregated over RTs to return only the distribution
#' of response and confidence ratings. Default: FALSE.
#'
#' @param stimulus numeric vector. Either 1, -1 or c(-1, 1) (default). Together with
#' condition represents the experimental situation. In a binary decision task the presented
#' stimulus belongs to one of two categories. In the default setting trials with
#' both categories presented are simulated but one can choose to simulate only trials with
#' the stimulus coming from one category (1 for the category that is associated with positive
#' drift in the decision process where "upper"/1 responses are considered correct and -1
#' correspondingly for negative drifts and "lower"/-1 correct decisions).
#'
#' @param seed numerical. Seeding for non-random data generation.
#' @param process_results logical. Whether the output simulations should contain the final
#' state of the decision (and visibility) process as additional column. Default is FALSE, meaning that
#' no additional columns for the final process states are returned.
#'
#' @return Depending on `gamma` and `agg_simus`.
#'
#' If `gamma` is `FALSE`, returns a `data.frame` with columns: `condition`,
#' `stimulus`, `response`, `correct`, `rt`, `conf` (the continuous confidence
#' measure) and `rating` (the discrete confidence rating), and `dec` and `vis`
#' (only if `process_results=TRUE`) for the final states of accumulators in the
#' simulation or
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
#' @details Simulation of response and decision times is done by simulating
#' normal variables in discretized steps until the lower or upper boundary
#' is met (or the maximal rt is reached). Afterwards, a confidence measure
#' is simulated according to the respective model.
#'
#' The confidence outputs are then binned according to the given thresholds.
#' The output of the fitting function \code{\link{fitRTConf}} with the respective model
#' fits the argument `paramDf` for simulation.
#' The Gamma coefficients are computed separately for correct/incorrect responses for the
#' correlation of confidence ratings with condition and rt and separately for conditions
#' for the correlation of accuracy and confidence. The
#' resulting data frames in the output thus have two columns. One for the grouping variable
#' and one for the Gamma coefficient.
#'
#' @note Different parameters for different conditions are only allowed for drift rate,
#' \code{v}, drift rate variability, \code{sv} and diffusion constant `s`.
#' All other parameters are used for all conditions.
#'
#' @references  Hellmann, S., Zehetleitner, M., & Rausch, M. (2023). Simultaneous modeling of choice, confidence and response time in visual perception. \emph{Psychological Review} 2023 Mar 13. doi: 10.1037/rev0000411. Epub ahead of print. PMID: 36913292.
#'
#'
#' @author Sebastian Hellmann.
#'
#' @name simulateWEV
#' @importFrom stats rnorm
#' @import dplyr
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @importFrom Rcpp evalCpp
#' @aliases simulate2DSD
#'
#' @examples
#' # Examples for "dynWEV" model (equivalent applicable
#' # for "2DSD" model (with less parameters))
#' # 1. Define some parameter set in a data.frame
#' paramDf <- data.frame(a=2.5,v1=0.1, v2=1, t0=0.1,z=0.55,
#'                       sz=0.3,sv=0.8, st0=0,  tau=3, w=0.1,
#'                       theta1=0.8, svis=0.5, sigvis=0.8)
#'
#' # 2. Simulate trials for both stimulus categories and all conditions (2)
#' simus <- simulateWEV(paramDf, model="dynWEV")
#' head(simus)
#' \donttest{
#'   library(ggplot2)
#'   simus <- simus[simus$response!=0,]
#'   simus$rating <- factor(simus$rating, labels=c("unsure", "sure"))
#'   ggplot(simus, aes(x=rt, group=interaction(correct, rating),
#'                     color=as.factor(correct), linetype=rating))+
#'     geom_density(size=1.2)+xlim(c(0,5))+
#'     facet_grid(rows=vars(condition), labeller = "label_both")
#' }
#'
#' # automatically aggregate simulation distribution
#' # to get only accuracy x confidence rating distribution for
#' # all conditions
#' agg_simus <- simulateWEV(paramDf, model="dynWEV", agg_simus = TRUE)
#' head(agg_simus)
#' \donttest{
#'   agg_simus$rating <- factor(agg_simus$rating, labels=c("unsure", "sure"))
#'   library(ggplot2)
#'   ggplot(agg_simus, aes(x=rating, group=correct, fill=as.factor(correct), y=p))+
#'     geom_bar(stat="identity", position="dodge")+
#'     facet_grid(cols=vars(condition), labeller = "label_both")
#' }
#' \donttest{
#'   # Compute Gamma correlation coefficients between
#'   # confidence and other behavioral measures
#'   # output will be a list
#'   simu_list <- simulateWEV(paramDf,n = 400, model="dynWEV", gamma=TRUE)
#'   simu_list
#' }

## When given vectorised parameters, n is the number of replicates for each parameter set
#' @rdname simulateWEV
#' @export
simulateWEV <- function (paramDf, n=1e+4,  model = "dynWEV", simult_conf = FALSE, gamma = FALSE, agg_simus=FALSE,
                         stimulus = c(-1,1), delta=0.01, maxrt=15, seed=NULL,
                         process_results = FALSE)
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
  if (!("lambda" %in% names(paramDf))) {
    if (model %in% c("dynaViTE", "2DSDT")) warning("No lambda specified in paramDf. lambda=0 used")
    lambda <- 0
  } else {
    lambda <- paramDf$lambda
  }
  if (model %in% c("WEVmu", "dynWEV")) model <- "dynaViTE"
  if (grepl("2DSD", model)) model <- "2DSD"
  if ("model" %in% names(paramDf)) paramDf$model <- NULL

  if (!(model %in% c("dynaViTE", "2DSD"))) stop("Only models dynaViTE, dynWEV (alias: WEVmu), and 2DSD/2DSDT are allowed!")

  if (!(all(stimulus %in% c(-1, 1)))) {
    stop(paste("Not accepted value for stimulus: ", paste(stimulus, collapse=", "),". Must be either 1, -1 or c(-1,1).", sep=""))
  }


  ## recover parameters from paramDf
  a <- paramDf$a
  z <- paramDf$z
  sz <- paramDf$sz
  t0 <- paramDf$t0  # recalc_t0 (see e.g. dWEV)
  st0 <- paramDf$st0
  tau = paramDf$tau
  if ("d" %in% names(paramDf)) {
    d <- paramDf$d
  } else {
    d <- 0
  }

  nConds <- length(grep(pattern = "^v[0-9]", names(paramDf), value = T))
  if (nConds > 0 ) {
    V <- c(t(paramDf[,paste("v",1:(nConds), sep = "")]))
  } else {
    V <- paramDf$v
    nConds <- 1
  }
  vary_sv <-   length(grep(pattern = "^sv[0-9]", names(paramDf), value = T))>1
  if (vary_sv){
    SV <- c(t((paramDf[,paste("sv",1:(nConds), sep = "")])))
  } else {
    SV <- rep(paramDf$sv, nConds)
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
  symmetric_confidence_thresholds <- length(grep(pattern = "thetaUpper", names(paramDf), value = T))<1
  if (symmetric_confidence_thresholds) {
    nRatings <- length(grep(pattern = "^theta[0-9]", names(paramDf)))+1
  } else {
    nRatings <- length(grep(pattern = "^thetaUpper[0-9]", names(paramDf)))+1
  }

  if (model =="2DSD") {
    w = -1
    svis = -1
    sigvis = -1
    muvis = rep(-1, nConds)
  } else {
    w = paramDf$w
    svis = paramDf$svis
    sigvis = paramDf$sigvis
    if ("muvis" %in% names(paramDf)) {
      muvis <- rep(paramDf$muvis, nConds)
    } else {
      muvis <- abs(V)
    }
  }

  df <- expand.grid(condition = 1:nConds, stimulus=stimulus)
  ## Produce process outcomes and compute confidence measure
  simus <- data.frame()
  for ( i in 1:nrow(df)) {
    s = as.numeric(S[df[i,]$condition])
    temp <- as.data.frame(r_WEV(n=n,
                                params=c(a/s,
                                         as.numeric(V[df[i,]$condition]*df[i,]$stimulus)/s,
                                         t0+st0/2, d, sz, as.numeric(SV[df[i,]$condition])/s,
                                         st0, z, tau,
                                         lambda,
                                         # when dynaViTE, use w, sigvis, svis, else do not supply the parameters
                                         rep(c(w, as.numeric(muvis[df[i,]$condition])/s,
                                               sigvis/s, svis/s),
                                             as.numeric(model=="dynaViTE"))),
                                delta = delta, maxT =maxrt, stop_on_error=TRUE))
    names(temp) <- c("rt", "response", "conf", "dec", "vis", "mu")
    temp$conf <- temp$conf * s
    temp$dec  <- temp$dec * s
    temp$vis  <- temp$vis * s

    simus <- rbind(simus,
                   cbind(condition=df[i, "condition"], stimulus=df[i, "stimulus"], temp))
  }
  if (!process_results) simus <- select(simus, -c("dec", "vis", "mu"))

  if (symmetric_confidence_thresholds) {
    thetas_upper <- c(-Inf, t(paramDf[,paste("theta",1:(nRatings-1), sep = "")]), Inf)
    thetas_lower <- c(-Inf, t(paramDf[,paste("theta",1:(nRatings-1), sep = "")]), Inf)
  } else {
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
  simus$rating[simus$response==-1] <- as.numeric(as.character(cut(simus$conf[simus$response==-1],
                                                                  breaks=thetas_lower, labels = levels_lower)))

  simus$correct <- as.numeric(simus$response == simus$stimulus)
  if (simult_conf) {
    simus$rt <- simus$rt + tau
  }

  if (gamma==TRUE) {
    gamma_condition <- simus %>% group_by(.data$correct) %>%
      summarise(Gamma=c(t(Hmisc::rcorr.cens(.data$rating,S=.data$condition, outx=TRUE)))[2],.groups="drop")
    gamma_rt <- simus %>% group_by(.data$correct) %>%
      summarise(Gamma=c(t(Hmisc::rcorr.cens(.data$rating,S=.data$rt, outx=TRUE)))[2],.groups="drop")
    gamma_correct <- simus %>% group_by(.data$condition) %>%
      summarise(Gamma=c(t(Hmisc::rcorr.cens(.data$rating,S=.data$correct, outx=TRUE)))[2],.groups="drop")
    gamma_rt_bycondition <- simus %>% group_by(.data$condition) %>%
      summarise(Gamma=c(t(Hmisc::rcorr.cens(.data$rating,S=.data$rt, outx=TRUE)))[2],.groups="drop")
    gamma_rt_byconditionbycorrect <- simus %>% group_by(.data$condition, .data$correct) %>%
      summarise(Gamma=c(t(Hmisc::rcorr.cens(.data$rating,S=.data$rt, outx=TRUE)))[2],.groups="drop")
  }
  if (agg_simus) {
    simus <- simus %>% group_by(.data$rating, .data$correct, .data$condition) %>%
      summarise(p = n()/(2*n),.groups="drop") %>%
      full_join(y=expand.grid(rating=1:nRatings, condition=1:nConds,
                              correct=c(0,1)), by=join_by("rating", "condition", "correct")) %>%
      mutate(p = ifelse(is.na(.data$p), 0, .data$p))
  } else {
    if (process_results) {
      simus <- simus[c("condition", "stimulus", "response", "correct", "rt", "dec", "vis", "mu", "conf", "rating")]
    } else {
      simus <- simus[c("condition", "stimulus", "response", "correct", "rt", "conf", "rating")]
    }
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
