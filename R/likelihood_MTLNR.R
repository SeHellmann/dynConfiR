#' Log-Likelihood functions for the mutliple-treshold log-normal race model
#'
#' Computes the Log-likelihood for given data and parameters in the MTLNR.
#' It is a wrapped version of the respective density \code{\link{dMTLNR}},
#' where one can find more information about the parameters. .
#' The function is mainly used inside \code{\link{fitRTConf}} for the MTLNR but exported
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
#'  If only one condition is used \code{v} will be accepted as well as \code{v1}.
#'
#'  \strong{stimulus, response}. stimulus and response should always
#'  be given in numerical format with values 1 and 2.
#'  Stimulus determines which of two accumulators has positive drift.
#'  The other has negative drift with the same absolute
#'  value. Response gives the index of the accumulator that reaches the
#'  boundary first.
#'
#' @references  Reynolds, A., Kvam, P. D., Osth, A. F., & Heathcote, A. (2020). Correlated racing evidence accumulator models. \emph{Journal of Mathematical Psychology, 96}, 102331. doi: doi: 10.1016/j.jmp.2020.102331
#'
#' @author Sebastian Hellmann.
#'
#' @name LogLikMTLNR
#' @importFrom dplyr case_when mutate rename
#' @importFrom magrittr %>%
#' @importFrom rlang .data
# @importFrom pracma integral
#' @aliases LLMTLNR
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
#' data <- rMTLNR(200,
#'                mu_v1 = as.numeric(stimulus==1)*discriminability*0.5,
#'                mu_v2 = as.numeric(stimulus==2)*discriminability*0.5,
#'                mu_d1=1, mu_d2=1, t0=0.1)
#' # discretize confidence ratings (only 2 steps: unsure vs. sure)
#' data$rating <- as.numeric(cut(data$conf, breaks = c(0, 3, Inf), include.lowest = TRUE))
#' data$stimulus <- stimulus
#' data$discriminability <- discriminability
#' data <- data[,-c(3,4)] # drop Tdec and conf measure (unobservable variable)
#' head(data)
#'
#' # 2. Define some parameter set in a data.frame
#' paramDf <- data.frame(v1=0.5, v2=1.0, t0=0.1, st0=0,
#'                       mu_d1=1, mu_d2=1,
#'                       s1=0.5, s2=0.5,
#'                       rho=0.2, theta1=2.5)
#'
#' # 3. Compute log likelihood for parameter and data
#' LogLikMTLNR(data, paramDf, condition="discriminability")
#'

#' @rdname LogLikMTLNR
#' @export
LogLikMTLNR <- function(data, paramDf,
                     precision=6, data_names = list(), ...) {
  #### Check data formatting ####
  tryCatch(data <- rename(data, ...),
           error = function(e) stop(paste0("Error renaming data columns. Probably a column name does not exist, or we tried to overwrite an already existing column.\nCheck whether an argument was misspelled and data name pairs are given in the form expected_name = true_name.\nUsed input for renaming columns:\n", paste(names(list(...)), list(...), sep="=", collapse = ", "))))


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


  ## Extract confidence thresholds
  if (symmetric_confidence_thresholds) {
    thetas_1 <- c(0, t(paramDf[,paste("theta",1:(nRatings-1), sep = "")]), 1e+128)
    thetas_2 <- c(0, t(paramDf[,paste("theta",1:(nRatings-1), sep = "")]), 1e+128)
  } else {
    thetas_1 <- c(0, t(paramDf[,paste("thetaUpper",1:(nRatings-1), sep = "")]), 1e+128)
    thetas_2 <- c(0, t(paramDf[,paste("thetaLower",1:(nRatings-1), sep="")]), 1e+128)
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
  ## Compute the row-wise mean pars for the accumulation rate:
  # For the stimulus accumulator, the accumulation rate is equal to the intensity
  # parameter v, for the other accumulator the mean accumulation rate is always 1
  # (as it scales the V parameters and mean boundary distances)
  data <-data %>% mutate(mu1 = V[.data$condition]*as.numeric(.data$stimulus==1),
                         mu2 = V[.data$condition]*as.numeric(.data$stimulus==2),
                         th1 = case_when(.data$response==1 ~ thetas_1[.data$rating],
                                         .data$response==2 ~ thetas_2[(.data$rating)]),
                         th2 = case_when(.data$response==1 ~ thetas_1[(.data$rating+1)],
                                         .data$response==2 ~ thetas_2[(.data$rating+1)]))

  if ("s_v1" %in% names(paramDf)) paramDf$s1 <- paramDf$s_v1
  if ("s_v2" %in% names(paramDf)) paramDf$s2 <- paramDf$s_v2
  if ("rho_v" %in% names(paramDf)) paramDf$rho <- paramDf$rho_v
  probs <- dMTLNR(
    data$rt, data$response, data$th1, data$th2,
     mu_v1=data$mu1, mu_v2=data$mu2,
     s_v1=paramDf$s1, s_v2=paramDf$s2,
     rho_v=paramDf$rho,
     mu_d1=paramDf$mu_d1, mu_d2=paramDf$mu_d2,
     s_d1=0, s_d2=0,
     rho_d=0,
     t0=paramDf$t0, st0=paramDf$st0,
     precision=precision)

  ## Produce output as log-Likelihood
  # if (any(is.na(probs))) return(-1e12)
  # if (any(probs<=0)) {
  #   return(-1e12)
  # }
  probs[probs==0 | is.na(probs)] <- .Machine$double.xmin

  if ("n" %in% names(data)) {
    logl <- sum(log(probs)*data$n)
  } else {
    logl <- sum(log(probs))
  }
  logl
}

