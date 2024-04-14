#' Simulation of confidence ratings and RTs in sequential sampling confidence models
#'
#' Simulates the decision responses, reaction times and confidence measure
#' together with a discrete confidence judgment for the sequential sampling confidence model
#' specified by the argument \code{model}, given specific parameter constellations.
#' This function is a wrapper that calls the respective functions for diffusion based
#' models (dynaViTE and 2DSD: \code{\link{simulateWEV}}) and race models (IRM, PCRM,
#' IRMt, and PCRMt: \code{\link{simulateRM}}. It also computes the Gamma rank correlation
#' between the confidence ratings and
#' condition (task difficulty), reaction times and accuracy in the simulated output.
#'
#' @param paramDf a list or dataframe with one row with the required parameters.
#' @param n integer. The number of samples (per condition and stimulus direction) generated.
#' Total number of samples is \code{n*nConditions*length(stimulus)}.
#' @param model character scalar. One of "dynaViTE", "dynWEV", "2DSD", "2DSDT",
#' "IRM", "PCRM", "IRMt", or "PCRMt". Could also be passed as a column in the paramDf argument.
#' @param gamma logical. If TRUE, the gamma correlation between confidence ratings, rt and accuracy is
#' computed.
#' @param agg_simus logical. Simulation is done on a trial basis with RTs outcome. If TRUE,
#' the simulations will be aggregated over RTs to return only the distribution of response and
#' confidence ratings. Default: FALSE.
#' @param  simult_conf logical. Whether in the experiment confidence was reported simultaneously
#' with the decision. If that is the case decision and confidence judgment are assumed to have happened
#' subsequent before the response. Therefore `tau` is included in the response time. If the decision was
#' reported before the confidence report, `simul_conf` should be `FALSE`.
#' @param stimulus numeric vector. Either 1, 2 or c(1, 2) (default).
#' Together with condition represents the experimental situation. In a 2AFC task the presented
#' stimulus belongs to one of two categories. In the default setting trials with
#' both categories presented are simulated but one can choose to simulate only trials with the
#' stimulus coming from one category.
#' @param delta numerical. Size of steps for the discretized simulation.
#' @param maxrt numerical. Maximum reaction time to be simulated. Default: 15.
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
#' @details The output of the fitting function \code{\link{fitRTConf}} with the respective model
#' fits the argument `paramDf` for simulation. The function calls the respective simulation
#' function for diffusion based models, i.e. dynaViTE and 2DSD (\code{\link{simulateWEV}}) or race models,
#' i.e. IRM(t) and PCRM(t), (\code{\link{simulateRM}}). See there for more information.
#'
#' \strong{Simulation Method:} The simulation is done by simulating normal variables
#' in discretized steps until
#' the processes reach the boundary. If no boundary is met within the maximum time,
#' response is set to 0.
#'
#' \strong{Gamma correlations:} The Gamma coefficients are computed separately for
#' correct/incorrect responses for the correlation of confidence ratings with condition and rt
#' and separately for conditions for the correlation of accuracy and confidence. The resulting
#' data frames in the output thus have two columns. One for the grouping variable and one for the
#' Gamma coefficient.
#'
#'
#' @author Sebastian Hellmann.
#'
#' @name simulateRTConf
#' @importFrom stats runif
# @importFrom pracma integral
#' @aliases rConfModel  simulateConfModel
#'
#' @examples
#'
#' # The function is particularly useful, when having a collection
#' # of parameter sets for different models (e.g. output by fitRTConfModels for
#' # more than one model).
#' library(dplyr)
#' # 1. Generate only one parameter set but for two different models
#' paramDf1 <- data.frame(model="dynWEV", a=1.5,v1=0.2, v2=1, t0=0.1,z=0.52,
#'                       sz=0.3,sv=0.4, st0=0,  tau=3, w=0.5,
#'                       theta1=1, svis=0.5, sigvis=0.8)
#' paramDf2 <- data.frame(model="PCRMt", a=2,b=2, v1=0.5, v2=1, t0=0.1,st0=0,
#'                       wx=0.6, wint=0.2, wrt=0.2, theta1=4)
#' paramDf <- full_join(paramDf1, paramDf2)
#' paramDf  # each model parameters sets hat its relevant parameters
#' # Split paramDf by model (maybe also other columns) and simulate data
#' simus <- paramDf |> group_by(model) |>
#'  reframe(simulateRTConf(cbind(cur_group(), pick(everything())), n=200, simult_conf = TRUE))
#' head(simus)
#'


## When given vectorised parameters, n is the number of replicates for each parameter set
#' @rdname simulateRTConf
#' @export
simulateRTConf <- function (paramDf, n=1e+4,  model = NULL,
                 gamma = FALSE, agg_simus=FALSE, simult_conf = FALSE,
                 stimulus = c(1,2), delta=0.01, maxrt=15, seed=NULL)
{
  gc(verbose = FALSE, full=FALSE)
  if (nrow(paramDf)>1) stop("paramDf must have one row.")
  paramDf <- paramDf[,c(!is.na(paramDf))]
  if (is.null(model) && ("model" %in% names(paramDf))) model <- paramDf$model
  if ((model %in% c("dynaViTE", "dynWEV", "WEVmu", "2DSD", "2DSDT")) && identical(stimulus, c(1,2))) stimulus <- c(-1,1)
  if (grepl("RM", model)) {
    res <- simulateRM(paramDf, n, model, FALSE, gamma, agg_simus, stimulus, delta, maxrt, seed)
    if (!agg_simus ) {
      if (gamma) {
        res$simus <- res$simus[c("condition", "stimulus", "response", "correct", "rt", "conf", "rating")]
      } else {
        res <- res[c("condition", "stimulus", "response", "correct", "rt", "conf", "rating")]
      }
    }
  } else if (model %in% c("dynaViTE", "dynWEV", "2DSD", "2DSDT")) {
    res <- simulateWEV(paramDf, n, model, simult_conf, gamma, agg_simus, stimulus, delta=delta, maxrt=maxrt, seed=seed)
  } else {stop("model not known.")}
  return(res)
}
