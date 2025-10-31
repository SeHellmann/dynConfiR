#' Simulation of confidence ratings and RTs in the
#' Multiple-threshold Log-normal Race Model
#'
#' Simulates the decision responses and reaction times together with a
#' discrete confidence judgment in the MTLNR (Reynolds et al., 2020),
#' given specific parameter constellations.
#' See \code{\link{dMTLNR}} for more information about parameters.
#' Also computes the Gamma rank correlation between the confidence ratings and condition
#' (task difficulty), reaction times and accuracy in the simulated output.
#' Basically, this function is a wrapper for \code{\link{rMTLNR}}
#' for application in confidence experiments with manipulation of specific parameters.
#'
#' @param paramDf a list or data frame with one row. Column names should
#' match the following names (see \link{dMTLNR}):
#' For different stimulus quality/mean
#' drift rates, names should be `v1`, `v2`, `v3`,.... (corresponding to the mean
#' parameter for the accumulation rate for the stimulus-corresponding accumulator,
#' therefore `mu_v1` and `mu_v2` are not used in this context but
#' replaced by the parameter `v`); `mu_d1` and `mu_d2` correspond to the mean
#' parameters for boundary distance of the two accumulators;
#' `s1` and `s2` correspond to the variance parameters of the first and
#' second boundary hitting time;
#' `rho` corresponds to the correlation of boundary hitting times.
#' Note that `s_v1`,`s_v2`,`rho_v`,`s_d1`,`s_d2`, and `rho_d` are not used in this
#' context, although the accumulation rate-related parameters can be used to replace
#' the above-mentioned variance parameters.
#' Additionally, the confidence thresholds should be given by names with
#' `thetaUpper1`, `thetaUpper2`,..., `thetaLower1`,... or,
#' for symmetric thresholds only by `theta1`, `theta2`,.... (see Details for the correspondence to the data)
#' @param n integer. The number of samples (per condition and stimulus direction) generated.
#' Total number of samples is \code{n*nConditions*length(stimulus)}.
#' @param gamma logical. If TRUE, the gamma correlation between confidence ratings, rt
#' and accuracy is computed.
#'
#' @param agg_simus logical. Simulation is done on a trial basis with RTs outcome.
#' If TRUE, the simulations will be aggregated over RTs to return only the distribution
#' of response and confidence ratings. Default: FALSE.
#'
#' @param stimulus numeric vector. Either 1, 2 or c(1,2) (default). Together with
#' condition represents the experimental situation. In a binary decision task the presented
#' stimulus belongs to one of two categories. In the default setting trials with
#' both categories presented are simulated but one can choose to simulate only trials with
#' the stimulus coming from one category (1 for the category that is associated with positive
#' drift in the decision process where 1 responses are considered correct and 2
#' correspondingly for negative drifts and 2 correct decisions).
#'
#' @param seed numerical. Seeding for non-random data generation.
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
#' @details Simulation is done by simulating normally distributed logarithms of
#' boundary crossing times for both accumulators based on the MTLNR model.
#' The smaller time determines decision time and response (i.e. the winning
#' accumulator). The confidence variable is computed based on the log-ratio of the
#' loosing boundary crossing time over the winning boundary crossing time.
#'
#' The confidence values are then binned according to the given thresholds.
#' The output of the fitting function \code{\link{fitRTConf}} with the respective model
#' fits the argument `paramDf` for simulation.
#' The Gamma coefficients are computed separately for correct/incorrect responses for the
#' correlation of confidence ratings with condition and rt and separately for conditions
#' for the correlation of accuracy and confidence. The
#' resulting data frames in the output thus have two columns. One for the grouping variable
#' and one for the Gamma coefficient.
#'
#' @note Different parameters for different conditions are only allowed for drift rate,
#' \code{v}.
#' All other parameters are used for all conditions.
#'
#' @references  Reynolds, A., Kvam, P. D., Osth, A. F., & Heathcote, A. (2020). Correlated racing evidence accumulator models. \emph{Journal of Mathematical Psychology, 96}, 102331. doi: doi: 10.1016/j.jmp.2020.102331
#'
#'
#' @author Sebastian Hellmann.
#'
#' @name simulateMTLNR
#' @import dplyr
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @importFrom Rcpp evalCpp
#'
#' @examples
#' # 1. Define some parameter set in a data.frame
#' paramDf <- data.frame(v1=0.5, v2=1.0, t0=0.1, st0=0,
#'                       mu_d1=1, mu_d2=1,
#'                       s1=0.5, s2=0.5, rho=0.2,
#'                       theta1=0.8, theta2=1.5)
#'
#' # 2. Simulate trials for both stimulus categories and all conditions (2)
#' simus <- simulateMTLNR(paramDf)
#' head(simus)
#' \donttest{
#'   library(ggplot2)
#'   simus <- simus[simus$response != 0, ]
#'   simus$rating <- factor(simus$rating, labels = c("unsure", "medium", "sure"))
#'   ggplot(simus, aes(x = rt, group = interaction(correct, rating),
#'                     color = as.factor(correct), linetype = rating)) +
#'     geom_density(linewidth = 1.2) + xlim(c(0, 5)) +
#'     facet_grid(rows = vars(condition), labeller = "label_both")
#' }
#'
#' # automatically aggregate simulation distribution
#' # to get only accuracy x confidence rating distribution for
#' # all conditions
#' agg_simus <- simulateMTLNR(paramDf, agg_simus = TRUE)
#' head(agg_simus)
#' \donttest{
#'   agg_simus$rating <- factor(agg_simus$rating, labels = c("unsure", "medium", "sure"))
#'   library(ggplot2)
#'   ggplot(agg_simus, aes(x = rating, group = correct, fill = as.factor(correct), y = p)) +
#'     geom_bar(stat = "identity", position = "dodge") +
#'     facet_grid(cols = vars(condition), labeller = "label_both")
#' }
#' \donttest{
#'   # Compute Gamma correlation coefficients between
#'   # confidence and other behavioral measures
#'   # output will be a list
#'   simu_list <- simulateMTLNR(paramDf, n = 400, gamma = TRUE)
#'   simu_list
#' }
#'
#' # Example with asymmetric confidence thresholds
#' paramDf_asym <- data.frame(v1=0.5, v2=1.0, t0=0.1, st0=0,
#'                           mu_d1=1, mu_d2=1,
#'                           s1=0.5, s2=0.5, rho=0.2,
#'                           thetaLower1=0.5, thetaLower2=1.2,
#'                           thetaUpper1=0.7, thetaUpper2=1.8)
#'
#' simus_asym <- simulateMTLNR(paramDf_asym, n = 1000)
#' head(simus_asym)
#'
#' # Example with multiple conditions
#' paramDf_multi <- data.frame(v1=0.3, v2=0.6, v3=1.2, t0=0.1, st0=0,
#'                            mu_d1=1, mu_d2=1,
#'                            s1=0.5, s2=0.5, rho=0.2,
#'                            theta1=0.8, theta2=1.5)
#'
#' simus_multi <- simulateMTLNR(paramDf_multi, n = 1000)
#' table(simus_multi$condition, simus_multi$correct)
#'

## When given vectorised parameters, n is the number of replicates for each parameter set
#' @rdname simulateMTLNR
#' @export
simulateMTLNR <- function (paramDf, n=1e+4, gamma = FALSE, agg_simus=FALSE,
                         stimulus = c(1,2), seed=NULL)
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

  if (!(all(stimulus %in% c(1, 2)))) {
    stop(paste("Not accepted value for stimulus: ", paste(stimulus, collapse=", "),". Must be either 1, 2, or c(1, 2).", sep=""))
  }

  ## recover parameters from paramDf

  nConds <- length(grep(pattern = "^v[0-9]", names(paramDf), value = T))
  if (nConds > 0 ) {
    V <- c(t(paramDf[,paste("v",1:(nConds), sep = "")]))
  } else {
    V <- paramDf$v
    nConds <- 1
  }
  symmetric_confidence_thresholds <- length(grep(pattern = "thetaUpper", names(paramDf), value = T))<1
  if (symmetric_confidence_thresholds) {
    nRatings <- length(grep(pattern = "^theta[0-9]", names(paramDf)))+1
  } else {
    nRatings <- length(grep(pattern = "^thetaUpper[0-9]", names(paramDf)))+1
  }

  if (symmetric_confidence_thresholds) {
    if (nRatings==1) {
      thresholds <- NULL
    } else {
      thresholds <- rep(c(t(paramDf[,paste("theta",1:(nRatings-1), sep = "")])), 2)
    }
  } else {
    thresholds <- c(t(paramDf[,paste("thetaLower",1:(nRatings-1), sep = "")]),
                    t(paramDf[,paste("thetaUpper",1:(nRatings-1), sep="")]))
  }

  if ("s_v1" %in% names(paramDf)) paramDf$s1 <- paramDf$s_v1
  if ("s_v2" %in% names(paramDf)) paramDf$s2 <- paramDf$s_v2
  if ("rho_v" %in% names(paramDf)) paramDf$rho <- paramDf$rho_v

  df <- expand.grid(condition = 1:nConds, stimulus=stimulus)
  ## Produce process outcomes and compute confidence measure
  simus <- data.frame()
  for ( i in 1:nrow(df)) {
    temp <- rMTLNR(n=n,thresholds=thresholds,
                 mu_v1=V[df[i,]$condition]*as.numeric(df[i,]$stimulus==1),
                 mu_v2=V[df[i,]$condition]*as.numeric(df[i,]$stimulus==2),
                 s_v1=paramDf$s1, s_v2=paramDf$s2,
                 rho_v=paramDf$rho,
                 mu_d1=paramDf$mu_d1, mu_d2=paramDf$mu_d2,
                 s_d1=0, s_d2=0, rho_d=0,
                 t0  = paramDf$t0, st0 = paramDf$st0)[,c("rt", "response", "conf", "rating")]
    simus <- rbind(simus,
                   cbind(condition=df[i, "condition"], stimulus=df[i, "stimulus"], temp))
  }
  simus$correct <- as.numeric(simus$response == simus$stimulus)

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
    simus <- simus[c("condition", "stimulus", "response", "correct", "rt", "conf", "rating")]
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
