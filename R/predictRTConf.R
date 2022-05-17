#' Prediction of confidence rating and response time distribution for sequential sampling confidence models
#'
#' \code{predictConf} predicts the categorical response distribution of
#' decision and confidence ratings, \code{predictRT} computes the predicted
#' RT distribution (density) for the sequential sampling confidence model
#' specified by the argument \code{model}, given specific parameter constellations.
#' This function calls the respective functions for diffusion based
#' models (dynWEV and 2DSD: \code{\link{predictWEV}}) and race models (IRM, PCRM,
#' IRMt, and PCRMt: \code{\link{predictRM}}).
#'
#' @param paramDf a list or dataframe with one row. Column names should match the
#' names of the respective model parameters. For different stimulus
#' quality/mean drift rates, names should be v1, v2, v3,.... Different s parameters
#' are possible with s1, s2, s3... with equally many steps as for drift rates (same
#' for sv parameter in dynWEV and 2DSD).
#' Additionally, the confidence thresholds should be given by names with
#' thetaUpper1, thetaUpper2,..., thetaLower1,... or,
#' for symmetric thresholds only by theta1, theta2,....
#' @param model character scalar. One of "2DSD", "dynWEV", "IRM", "PCRM", "IRMt", or "PCRMt".
#' @param maxrt numeric. The maximum RT for the
#' integration/density computation. Default: 15 (for \code{predictConf} (integration)),
#' 9 (for \code{predictRT}).
#' @param subdivisions \code{integer} (default: 100).
#' For \code{predictConf} it is used as argument for the inner integral routine.
#' For \code{predictRT} it is the number of points for which the density is computed.
#' @param minrt numeric or NULL(default). The minimum rt for the density computation.
#' @param  simult_conf logical, only relevant for dynWEV and 2DSD. Whether in the experiment
#' confidence was reported simultaneously with the decision, as then decision and confidence
#' judgment are assumed to have happened subsequent before response and computations are
#' different, when there is an observable interjudgment time (then `simult_conf` should be FALSE).
#' @param scaled logical. For \code{predictRT}. Whether the computed density
#' should be scaled to integrate to one (additional column densscaled). Otherwise the output
#' is a defective density (i.e. its integral is equal to the probability of a response and
#' not 1). If TRUE, the argument `DistConf` should be given, if available. Default: FALSE.
#' @param DistConf NULL or data.frame. For \code{predictRT}. A data.frame or matrix
#' with column names, giving the distribution of response and rating choices for
#' different conditions and stimulus categories in the form of the output of
#' \code{predictConf}. It is only necessary, if `scaled=TRUE`, because these
#' probabilities are used for scaling. If `scaled=TRUE` and `DistConf=NULL`, it will be
#' computed with the function \code{predictConf}, which takes some time and the function will
#' throw a message. Default: NULL
#' @param stop.on.error logical. Argument directly passed on to integrate. Default is FALSE,
#' since the densities invoked may lead to slow convergence of the integrals (which are still
#' quite accurate) which causes R to throw an error.
#' @param .progress logical. If TRUE (default) a progress bar is drawn to the console.
#'
#' @return \code{predictConf} gives a data frame/tibble with columns: condition, stimulus,
#' response, rating, correct, p, info, err. p is the predicted probability of a response
#' and rating, given the stimulus category and condition. Message and error refer to the
#' respective outputs of the integration routine used for computation.
#' \code{predictRT} returns a data frame/tibble with columns: condition, stimulus,
#' response, rating, correct, rt and dens (and densscaled, if `scaled=TRUE`).
#'
#'
#' @details The function \code{predictConf} consists merely of an integration of
#' the reaction time density of the given model, \code{{d*model*}}, over the response
#' time in a reasonable interval (0 to maxrt). The function \code{predictRT} wraps
#' these density functions to a parameter set input and a data.frame output.
#' For the argument \code{paramDf}, the output of the fitting function \code{\link{fitRTConf}}
#' with the respective model may be used.
#'
#' @note Different parameters for different conditions are only allowed for drift rate,
#' \code{v}, drift rate variability, \code{sv} (in dynWEV and 2DSD), and process variability
#' `s`. All other parameters are used for all conditions.
#'
#' @author Sebastian Hellmann.
#'
#' @name predictRTConf
#' @importFrom stats integrate
#' @import dplyr
#' @importFrom progress progress_bar
#' @importFrom magrittr %>%
#' @importFrom rlang .data
# @importFrom pracma integral
#' @aliases predictConf
#' @importFrom Rcpp evalCpp
#'

#' @rdname predictRTConf
#' @export
predictConf <- function(paramDf, model=NULL,
                           maxrt=15, subdivisions = 100L, simult_conf = FALSE, stop.on.error=FALSE,
                           .progress=TRUE){
  paramDf <- as.data.frame(paramDf)
  paramDf <- paramDf[,!is.na(paramDf)]
  if (is.null(model) && ("model" %in% names(paramDf))) model <- paramDf$model

  #### Check model argument
  if (grepl("RM", model)) {
    res <- predictRM_Conf(paramDf, model, FALSE, maxrt, subdivisions, stop.on.error, .progress)
  } else if (model %in% c("dynWEV", "2DSD")) {
    res <- predictWEV_Conf(paramDf, model, 3, maxrt, subdivisions,simult_conf, stop.on.error, .progress)
  } else { stop("model not known.")}
  return(res)
}



### Predict RT-distribution
#' @rdname predictRTConf
#' @export
predictRT <- function(paramDf, model=NULL,
                         maxrt=9, subdivisions = 100L,  minrt=NULL, simult_conf=FALSE,
                         scaled = FALSE, DistConf=NULL,
                         .progress = TRUE) {
  paramDf <- as.data.frame(paramDf)
  paramDf <- paramDf[,!is.na(paramDf)]
  if (is.null(model) && ("model" %in% names(paramDf))) model <- paramDf$model

  #### Check model argument
  if (grepl("RM", model)) {
    res <- predictRM_RT(paramDf, model, FALSE, maxrt, subdivisions, minrt, scaled, DistConf, .progress)
  } else if (model %in% c("dynWEV", "2DSD")) {
    res <- predictWEV_RT(paramDf, model, 3, maxrt, subdivisions, minrt, simult_conf, scaled, DistConf, .progress)
  } else { stop("model not known.")}
  return(res)
}