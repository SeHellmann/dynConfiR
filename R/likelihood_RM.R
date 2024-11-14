#' Log-Likelihood functions for the independent and partially anti-correlated race models of confidence
#'
#' Computes the Log-likelihood for given data and parameters in the IRM and PCRM with or without time-scaled
#' confidence measure. It is a wrapped version of the respective densities \code{\link{dIRM}} and \code{\link{dPCRM}},
#' where one can find more information about the parameters. It restricts the rates of accumulation to be the negative
#' of each other, though (a common assumption in perceptual decision tasks).
#' The function is mainly used inside \code{\link{fitRTConf}} for race models but exported
#' for individual usage in other contexts.
#'
#' @param data a dataframe where each row is one trial. Containing following variables:
#' \itemize{
#'   \item condition    (not necessary; convertible to integer (e.g. factor); for different levels of stimulus quality),
#'   \item rating            (convertible to integer (e.g. factor); discrete confidence judgments),
#'   \item rt                    (numeric; giving reaction times for decision task),
#'   \item stimulus     (values at least convertible to c(1,2), i.e. integer or factor; stimulus category (index of accumulator with higher drift))
#'   \item response     (values at least convertible to c(1,2); direction of decision; (index of accumulator reaching the boundary first))
#'
#' }
#' @param paramDf a list or data frame with one row. Column names should match the names of
#' \link{RaceModels} parameter names (only `mu1` and `mu2` are not used in this context but
#' replaced by the parameter `v`). For different stimulus quality/mean
#' drift rates, names should be `v1`, `v2`, `v3`,....
#' Different `s` parameters are possible with `s1`, `s2`, `s3`,... with equally many steps as for drift rates. Additionally, the confidence
#' thresholds should be given by names with `thetaUpper1`, `thetaUpper2`,..., `thetaLower1`,... or,
#' for symmetric thresholds only by `theta1`, `theta2`,.... (see Details for the correspondence to the data)
#' @param model character scalar. One of "IRM" or "PCRM". ("IRMt" and "PCRMt" will also be accepted. In that case,
#' time_scaled is set to TRUE.)
#' @param time_scaled logical. Whether the confidence measure should be scaled by 1/sqrt(rt). Default: TRUE.
#' @param precision numerical scalar. Precision of calculation for integration over t0.
#' @param data_names list. Possibility of giving alternative column names for the variables in the data. By default column names are identical to the
#' ones given in the data argument description.
#' @param ... Another possibility of giving alternative variable names in data frame (in the form \code{condition = "SOA"}).
#'
#' @return Numeric scalar. The summed Log-likelihood of the data given the parameters in the respective model. If one or more row-wise probabilities is <=0,
#' the function returns -1e+12.
#'
#' @details Note, that the requirements on the format of the columns for the likelihood functions are much stricter, than in \code{\link{fitRTConf}}.
#' This is because the function is very frequently called in the optimization routines of the fitting process and the preprocessing steps are
#' therefore included in the other function.
#'
#'  \strong{rating, condition}. If integer, values should range from 1
#'  to number of possible ratings/conditions. If factor,
#'  the number of levels should be equal to number of possible
#'  ratings/conditions. This should be consistent with the
#'  parameter vector. The confidence thresholds should be named as
#'  `thetaUpper1`, `thetaLower1`,.... (or `theta1`,... for symmetric
#'  thresholds), with the number of ratings -1 and the mean drift rates
#'  (and possibly the standard deviation in drift rates)
#'  should be denoted as `v1`, `v2`,...
#'  (and `s1`, `s2`,...) with the number equal to the number of conditions.
#'  If only one condition is used \code{v} will be accepted as well as \code{v1}.
#'
#'  \strong{stimulus, response}. stimulus and response should always
#'  be given in numerical format with values 1 and 2.
#'  Stimulus determines which of two accumulators has positive drift.
#'  The other has negative drift with the same absolute
#'  value. Response gives the index of the accumulator that reaches the
#'  boundary first.
#'
#'
#' @author Sebastian Hellmann.
#'
#' @name LogLikRM
#' @importFrom dplyr case_when mutate rename
#' @importFrom magrittr %>%
#' @importFrom rlang .data
# @importFrom pracma integral
#' @aliases LLRM LogLikIRM  LogLikPCRM
#' @importFrom Rcpp evalCpp
#'
#' @examples
#' # 1. Generate data from an artificial participants
#' # Get random index for accumulator with positive
#' # drift (i.e. stimulus category) and
#' # stimulus discriminability (two steps: hard, easy)
#' stimulus <- sample(c(1, 2), 200, replace=TRUE)
#' discriminability <- sample(c(1, 2), 200, replace=TRUE)
#' # generate data for participant 1
#' data <- rPCRM(200, mu1=ifelse(stimulus==1, 1, -1)*discriminability*0.5,
#'               mu2=ifelse(stimulus==1, -1, 1)*discriminability*0.5,
#'              a=2, b=1.8, t0=0.2, st0=0, wx=0.7, wint=0.3, wrt=0)
#' # discretize confidence ratings (only 2 steps: unsure vs. sure)
#' data$rating <- as.numeric(cut(data$conf, breaks = c(0, 3, Inf), include.lowest = TRUE))
#' data$stimulus <- stimulus
#' data$discriminability <- discriminability
#' data <- data[data$response!=0, ] # drop not finished decision processes
#' data <- data[,-c(3,4)] # drop xl and conf measure (unobservable variable)
#' head(data)
#'
#' # 2. Define some parameter set in a data.frame
#' paramDf <- data.frame(a=2,b=2, v1=0.5, v2=1, t0=0.1,st0=0,
#'                       wx=0.6, wint=0.2, wrt=0.2,
#'                       theta1=4)
#'
#' # 3. Compute log likelihood for parameter and data
#' LogLikRM(data, paramDf, model="PCRMt", condition="discriminability")
#' # same result
#' LogLikRM(data, paramDf, model="PCRM", time_scaled=TRUE,condition="discriminability")
#' # different
#' LogLikRM(data, paramDf, model="PCRM", condition="discriminability")
#'
#' # same parameters used for IRM model
#' LogLikRM(data, paramDf, model="IRMt", condition="discriminability")
#'

#' @rdname LogLikRM
#' @export
LogLikRM <- function(data, paramDf, model="IRM", time_scaled =FALSE,
                     precision=6, data_names = list(), ...) {
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

  #### Check data formatting ####
  # ToDo
  data <- rename(data, ...)


  #### Get information from paramDf ####
  nConds <- length(grep(pattern = "^v[0-9]", names(paramDf), value = T))
  symmetric_confidence_thresholds <- length(grep(pattern = "thetaUpper", names(paramDf), value = T))<1
  if (symmetric_confidence_thresholds) {
    nRatings <- length(grep(pattern = "^theta[0-9]", names(paramDf)))+1
  } else {
    nRatings <- length(grep(pattern = "^thetaUpper[0-9]", names(paramDf)))+1
  }
  if (nConds > 0 ) {
    V <- c(t(paramDf[,paste("v",1:(nConds), sep = "")]))
  } else {
    V <- paramDf$v
    nConds <- 1
  }


  ## Recover confidence thresholds
  if (symmetric_confidence_thresholds) {
    thetas_1 <- c(1e-32, t(paramDf[,paste("theta",1:(nRatings-1), sep = "")]), 1e+64)
    thetas_2 <- c(1e-32, t(paramDf[,paste("theta",1:(nRatings-1), sep = "")]), 1e+64)
  } else {
    thetas_1 <- c(1e-32, t(paramDf[,paste("thetaUpper",1:(nRatings-1), sep = "")]), 1e+64)
    thetas_2 <- c(1e-32, t(paramDf[,paste("thetaLower",1:(nRatings-1), sep="")]), 1e+64)
  }
  #### Check for column names given ####
  names_missing <- !(c("condition","response","stimulus","rating", "rt", "sbj", "correct") %in% names(data_names))
  data_names <- c(data_names,
                  setNames(as.list(c("condition","response","stimulus","rating", "rt", "sbj", "correct")[names_missing]),
                           c("condition","response","stimulus","rating", "rt", "sbj", "correct")[names_missing]))
  if (is.null(data[[data_names$condition]])) data[[data_names$condition]] <- 1

  if (!(all(data[[data_names$response]] %in% c(1,2)) & all(data[[data_names$stimulus]] %in% c(1,2)))) {
    data <- data %>% mutate(response = if_else(.data[[data_names$response]]==sort(unique(data[[data_names$response]]))[1],1,2),
                            stimulus = if_else(.data[[data_names$stimulus]]==sort(unique(data[[data_names$stimulus]]))[1],1,2))
  } else {
    data$response <- data[[data_names$response]]
    data$stimulus <- data[[data_names$stimulus]]
  }


  data$condition <- as.numeric(factor(data[[data_names$condition]],levels = sort(unique(data[[data_names$condition]]))))

  if (!is.numeric(data$rating)) {
    data$rating <- as.numeric(as.factor(data[[data_names$rating]]))
  }
  ## Compute the row-wise likelihood of observations
  data <-data %>% mutate(a = paramDf$a,
                         b = paramDf$b,
                         th1 = case_when(.data$response==1 ~ thetas_1[(.data$rating)],
                                         .data$response==2 ~ thetas_2[(.data$rating)]),
                         th2 = case_when(.data$response==1 ~ thetas_1[(.data$rating+1)],
                                         .data$response==2 ~ thetas_2[(.data$rating+1)]),
                         mu1 = V[.data$condition]*(-1)^(1+.data$stimulus),
                         mu2 = V[.data$condition]*(-1)^(.data$stimulus),
                         t0 = paramDf$t0,
                         st0 = paramDf$st0)
  if (time_scaled) {
    data$wx = paramDf$wx
    data$wrt = paramDf$wrt
    data$wint = paramDf$wint
  } else {
    data$wx = 1
    data$wrt = 0
    data$wint = 0
  }
  if ("s" %in% names(paramDf)) {
    data$s = paramDf$s
  } else {
    data$s = 1
  }
  if (model=="IRM") {
    probs <- dIRM(data$rt, data$response,data$mu1, data$mu2, data$a, data$b,
                  data$th1, data$th2, data$wx,  data$wrt,  data$wint,
                  data$t0, data$st0, data$s, time_scaled=time_scaled,
                  precision=precision)
  } else {
    probs <- dPCRM(data$rt, data$response,data$mu1, data$mu2, data$a, data$b,
                   data$th1, data$th2, data$wx,  data$wrt,  data$wint,
                   data$t0, data$st0, data$s, time_scaled=time_scaled,
                   precision=precision)
  }

  ## Produce output as log-Likelihood
  # if (any(is.na(probs))) return(-1e12)
  # if (any(probs<=0)) {
  #   return(-1e12)
  # }
  probs[probs==0] <- .Machine$double.xmin

  if ("n" %in% names(data)) {
    logl <- sum(log(probs)*data$n)
  } else {
    logl <- sum(log(probs))
  }
  logl
}

