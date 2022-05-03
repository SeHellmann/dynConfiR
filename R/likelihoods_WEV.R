#' Log-Likelihood functions for the dWEV and 2DSD models
#'
#' Computes the Log-likelihood for given data and parameters in the versions of the dWEV model (Hellmann & Rausch) and the 2DSD model (Pleskac & Busemeyer, 2010).
#' It is a wrapped version of the respective densities \code{\link{dWEV}} and \code{\link{d2DSD}}, where one can find more information about the parameters (\code{z} is always given relatively, in the likelihood).
#' The function is mainly used in \code{\link{fitRTConf}} but exported for individual usage in other contexts.
#'
#' @param data a dataframe where each row is one trial. Containing following variables:
#' \itemize{
#'   \item condition    (not necessary; convertable to integer (e.g. factor); for different levels of stimulus quality),
#'   \item rating            (convertable to integer (e.g. factor); discrete confidence judgments),
#'   \item rt                    (numeric; giving reaction times for decision task),
#'   \item stimulus     (values at least convertible to c(-1,1); stimulus category (direction of evidence accumulation))
#'   \item response     (characters in "upper", "lower" (or convertible to); direction of decision; correct answers are "lower" for stimulus=-1; and "upper" for stimulus=1),
#'
#' }
#' @param paramDf list or dataframe with one row. Names should match the names of dWEV and 2DSD model specific parameter names.
#' For different stimulus quality/mean drift rates, names should be v1, v2, v3,.... Different sv parameters are possible with sv1, sv2, sv3...
#' with equally many steps as for drift rates. Additionally, the confidence thresholds should be given by names with
#' thetaUpper1, thetaUpper2,..., thetaLower1,... or, for symmetric thresholds only by theta1, theta2,.... (see Details for the correspondance to the data)
#' @param model character scalar. One of "dynWEV" or "2DSD" for the model to fit.
#' @param  simult_conf logical. Whether in the experiment confidence was reported simultaneously
#' with the decision, as then decision and confidence judgment are assumed to have happened
#' subsequent before response and computations are different, when there is an observable
#' interjudgment time (then simult_conf should be FALSE).
#' @param data_names list. Possibility of giving alternative column names for the variables in the data. By default columnnames are identical to the
#' ones given in the data argument description.
#' @param precision numerical scalar. Precision of calculation for integration over z and t0.
#' @param stop_on_error logical. If TRUE an error in the function will be returned in case of invalid parameters. Otherwise, the output will be 0 without error.
#' @param ... Possibility of giving alternative variable names (unquoted) in data frame (in the form \code{condition = SOA}).
#'
#' @return Numeric scalar. The summed Log-likelihood of the data given the parameters in the respective model. If one or more row-wise probabilities is <=0,
#' the function returns -1e+12.
#'
#' @details Note, that the requirements on the format of the columns for the likelihood functions are much stricter, than in \code{\link{fitRTConf}}.
#' This is because the function is very frequently calls in the optimisation routines of the fitting process and the preprocessing steps are
#' therefore included in that function.
#'
#'  \strong{rating, condition}. If integer, values should range from 1 to number of possible ratings/conditions. If factor, the number of levels should be
#'  equal to number of possible ratings/conditions. This should be consistent with the parameter vector. The confidence thresholds should be named
#'  as thetaUpper1, thetaLower1,.... (or theta1,... for symmetric thresholds), with the number of ratings -1 and the mean drift rates (and possibly the
#'  standard deviation in drift rates) should be denoted as v1, v2,... (and sv1, sv2,...) with the number equal to the number of conditions. If only
#'  one condition is used \code{v} will be accepted as well as \code{v1}.
#'
#'  \strong{stimulus, response}. stimulus should always be given in numerical format with values -1 and 1.
#'  reponse should always be given as a character vector with "lower" and "upper" as values.
#'  This corresponds to the situation of Ratcliff's diffusion model (Ratcliff, 1978), where stimulus is the sign of the mean drift direction and
#'  the response is the "upper" or "lower" boundary that is first hit by the evidence accumulation. A correct decision is therefore "lower",
#'  if stimulus is -1, and "upper", if stimulus is 1.
#'

#' @references Pleskac, T. J., & Busemeyer, J. R. (2010). Two-Stage Dynamic Signal Detection: A Theory of Choice, Decision Time, and Confidence, \emph{Psychological Review}, 117(3), 864-901. doi:10.1037/a0019737
#'
#' Ratcliff, R. (1978). A theory of memory retrieval. \emph{Psychological Review}, 85(2), 59-108.
#'
#' Rausch, M., Hellmann, S., & Zehetleitner, M. (2018). Confidence in masked orientation judgments is informed by both evidence and visibility. \emph{Attention, Perception, & Psychophysics}, 80(1), 134–154.  doi: 10.3758/s13414-017-1431-5
#'
#'
#' @author Sebastian Hellmann.
#'
#' @name LogLikWEV
#' @importFrom dplyr case_when mutate rename
#' @importFrom magrittr %>%
#' @importFrom rlang .data
# @importFrom pracma integral
#' @aliases LLWEV LLdynWEV LL2DSD  LogLikWEV  LogLikdynWEV LogLikWEV2DSD
#' @importFrom Rcpp evalCpp
#'

#' @rdname LogLikWEV
#' @export
LogLikWEV <- function(data, paramDf, model="dynWEV", simult_conf = FALSE, precision=1e-5, stop_on_error = TRUE, data_names = list(), ...) {
  #### Check data formatting ####
  data <- rename(data, ...)

  #### Get information from paramDf ####
  nConds <- length(grep(pattern = "^v[0-9]", names(paramDf), value = T))
  symmetric_confidence_thresholds <- length(grep(pattern = "thetaUpper", names(paramDf), value = T))<1
  if (symmetric_confidence_thresholds) {
    nRatings <- length(grep(pattern = "^theta[0-9]", names(paramDf)))+1
  } else {
    nRatings <- length(grep(pattern = "^thetaUpper[0-9]", names(paramDf)))+1
  }
  vary_sv <-   length(grep(pattern = "^sv[0-9]", names(paramDf), value = T))>1

  if (nConds > 0 ) {
    V <- c(t(paramDf[,paste("v",1:(nConds), sep = "")]))
  } else {
    V <- paramDf$v
    nConds <- 1
  }
  if (vary_sv){
    SV <- c(t((paramDf[,paste("sv",1:(nConds), sep = "")])))
  } else {
    SV <- rep(paramDf$sv, nConds)
  }

  ## Recover confidence thresholds
  if (symmetric_confidence_thresholds) {
    thetas_upper <- c(-1e+32, t(paramDf[,paste("theta",1:(nRatings-1), sep = "")]), 1e+32)
    thetas_lower <- c(-1e+32, t(paramDf[,paste("theta",1:(nRatings-1), sep = "")]), 1e+32)
  } else {
    thetas_upper <- c(-1e+32, t(paramDf[,paste("thetaUpper",1:(nRatings-1), sep = "")]), 1e+32)
    thetas_lower <- c(-1e+32, t(paramDf[,paste("thetaLower",1:(nRatings-1), sep="")]), 1e+32)
  }
  # if (model=="2DSD") {    # For 2DSD the parametrisation for lower thetas is different (different confidence scale)
  #   thetas_lower <- c(-1e+32, rev(thetas_lower[2:(nRatings)]), 1e+32)
  #   if (symmetric_confidence_thresholds) {
  #     thetas_lower <- paramDf$a- rev(thetas_upper)
  #   }
  # }

  #### Check for column names given ####
  names_missing <- !(c("condition","response","stimulus","rating", "rt", "sbj", "correct") %in% names(data_names))
  data_names <- c(data_names,
                  setNames(as.list(c("condition","response","stimulus","rating", "rt", "sbj", "correct")[names_missing]),
                           c("condition","response","stimulus","rating", "rt", "sbj", "correct")[names_missing]))
  if (is.null(data[[data_names$condition]])) data[[data_names$condition]] <- 1

  data <- data %>% mutate(response = if_else(.data[[data_names$response]]==sort(unique(data[[data_names$response]]))[1],"lower","upper"),
                          stimulus = if_else(.data[[data_names$stimulus]]==sort(unique(data[[data_names$stimulus]]))[1],-1,1),
                          condition = as.numeric(factor(.data[[data_names$condition]],levels = sort(unique(data[[data_names$condition]])))))
  if (!is.numeric(data$rating)) {
    data <- data %>% mutate(rating = as.numeric(as.factor(.data[[data_names$rating]])))
  }
  ## Compute the row-wise likelihood of observations

  data <-data %>% mutate(vth1 = case_when(.data$response=="upper" ~ thetas_upper[.data$rating],
                                          .data$response=="lower" ~ thetas_lower[(.data$rating)]),
                         vth2 = case_when(.data$response=="upper" ~ thetas_upper[(.data$rating+1)],
                                          .data$response=="lower" ~ thetas_lower[(.data$rating+1)]),
                         M_drift = V[.data$condition]*.data$stimulus,
                         SV = SV[.data$condition])
  probs <- with(data, switch(which(model== c("dynWEV", "2DSD")),
                             dWEV(rt, vth1,vth2,
                                        response=response,
                                        tau=paramDf$tau, a=paramDf$a,
                                        v = M_drift,
                                        t0 = paramDf$t0, z = paramDf$z, sz = paramDf$sz, st0=paramDf$st0,
                                        sv = SV, w=paramDf$w, muvis=abs(M_drift), svis=paramDf$svis,
                                        sigvis=paramDf$sigvis,
                                        simult_conf = simult_conf, z_absolute = FALSE,
                                        precision = precision, stop_on_error = stop_on_error),
                             d2DSD(rt, vth1,vth2,
                                   response=response, tau=paramDf$tau, a=paramDf$a,
                                   v = M_drift,
                                   t0 = paramDf$t0, z = paramDf$z, sz = paramDf$sz, st0=paramDf$st0,
                                   sv = SV, simult_conf = simult_conf, z_absolute = FALSE,
                                   precision = precision, stop_on_error = stop_on_error)))



  ## Produce output as log-Likelihood
  if (any(is.na(probs))) return(-1e12)
  if (any(probs<=0)) {
    return(-1e12)
  }
  if ("n" %in% names(data)) {
    logl <- sum(log(probs)*data$n)
  } else {
    logl <- sum(log(probs))
  }
  logl
}
