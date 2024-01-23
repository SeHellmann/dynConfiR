#' Log-Likelihood functions for the dynWEV and 2DSD models of confidence
#'
#' Computes the Log-likelihood for given data and parameters in the
#' dynWEV model (Hellmann et al., 2023) and the 2DSD model
#' (Pleskac & Busemeyer, 2010). It is a wrapped version of the
#' respective densities \code{\link{dWEV}} and \code{\link{d2DSD}},
#' where one can find more information about the parameters
#' (\code{z} is always given relatively, in the likelihood).
#' The function is mainly used in \code{\link{fitRTConf}} but
#' exported for individual usage in other contexts.
#'
#' @param data a dataframe where each row is one trial. Containing following variables:
#' \itemize{
#'   \item condition    (not necessary; convertible to integer (e.g. factor); for different levels of stimulus quality),
#'   \item rating       (convertible to integer (e.g. factor); discrete confidence judgments),
#'   \item rt           (numeric; giving reaction times for decision task),
#'   \item stimulus     (values at least convertible to c(-1,1); stimulus category (direction of evidence accumulation))
#'   \item response     (characters in `"upper"`, `"lower"` (or convertible to); direction of decision; correct answers are "lower" for stimulus=-1; and "upper" for stimulus=1),
#'
#' }
#' @param paramDf list or data.frame with one row. Names should match the names of
#' \link{dynaViTE} and \link{2DSD} model specific parameter names. For different
#' stimulus quality/mean drift rates, names should be `v1`, `v2`, `v3`,....
#' Different `sv` and/or `s` parameters are possible with `sv1`, `sv2`, `sv3`... (`s1`, `s2`, `s3`,...
#' respectively) with equally many steps as for drift rates. Additionally, the confidence
#' thresholds should be given by names with `thetaUpper1`, `thetaUpper2`,..., `thetaLower1`,... or,
#' for symmetric thresholds only by `theta1`, `theta2`,...
#' (see Details for the correspondence to the data)
#' @param model character scalar. One of "dynWEV" or "2DSD" for the model to fit.
#' @param  simult_conf logical. Whether in the experiment confidence was reported simultaneously
#' with the decision, as then decision and confidence judgment are assumed to have happened
#' subsequent before response and computations are different, when there is an observable
#' interjudgment time (then `simult_conf` should be `FALSE`).
#' @param data_names list. Possibility of giving alternative column names for the variables
#' in the data. By default column names are identical to the
#' ones given in the data argument description.
#' @param precision numerical scalar. Precision of calculation for integration over z and t0.
#' @param stop_on_error logical. If TRUE an error in the function will be returned in case of
#' invalid parameters. Otherwise, the output will be 0 without error.
#' @param ... Possibility of giving alternative variable names in data frame
#' (in the form \code{condition = "SOA"}).
#'
#' @return Numeric scalar. The summed Log-likelihood of the data given the parameters in the respective model. If one or more row-wise probabilities is <=0,
#' the function returns -1e+12.
#'
#' @details Note, that the requirements on the format of the columns for the likelihood functions
#' are much stricter, than in \code{\link{fitRTConf}}.
#' This is because the function is very frequently calls in the optimization routines of the
#' fitting process and the preprocessing steps are
#' therefore included in that function.
#'
#'  \strong{rating, condition}. If integer, values should range from 1 to number of possible
#'  ratings/conditions. If a factor, the number of levels should be
#'  equal to number of possible ratings/conditions. This should be consistent with the
#'  parameter vector. The confidence thresholds should be named
#'  as `thetaUpper1`, `thetaLower1`,.... (or `theta1`,... for symmetric thresholds), with the
#'  number of ratings -1 and the mean drift rates (and possibly the
#'  standard deviation in drift rates) should be denoted as `v1`, `v2`,...
#'  (and `sv1`, `sv2`,.../`s1`, `s2`, ...) with the number equal to the number of conditions.
#'  If only one condition is used \code{v} will be accepted as well as \code{v1}.
#'
#'  \strong{stimulus, response}. stimulus should always be given in numerical format with values -1 and 1.
#'  response should always be given as a character vector with `"lower"` and `"upper"` as values.
#'  This corresponds to the situation of Ratcliff's diffusion model (Ratcliff, 1978), where stimulus is the sign of the mean drift direction and
#'  the response is the `"upper"` or `"lower"` boundary that is first hit by the evidence accumulation. A correct decision is therefore `"lower"`,
#'  if stimulus is -1, and `"upper"`, if stimulus is 1.
#'
#' @references  Hellmann, S., Zehetleitner, M., & Rausch, M. (2023). Simultaneous modeling of choice, confidence and response time in visual perception. \emph{Psychological Review} 2023 Mar 13. doi: 10.1037/rev0000411. Epub ahead of print. PMID: 36913292.
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
#' @examples
#' # 1. Generate data from an artificial participants
#' # Get random drift direction (i.e. stimulus category) and
#' # stimulus discriminability (two steps: hard, easy)
#' stimulus <- sample(c(-1, 1), 200, replace=TRUE)
#' discriminability <- sample(c(1, 2), 200, replace=TRUE)
#' # generate data for participant 1
#' data <- rWEV(200, a=2,v=stimulus*discriminability*0.5,
#'              t0=0.2,z=0.5, sz=0.1,sv=0.1, st0=0,  tau=4, s=1, w=0.3)
#' # discretize confidence ratings (only 2 steps: unsure vs. sure)
#' data$rating <- as.numeric(cut(data$conf, breaks = c(-Inf, 1, Inf), include.lowest = TRUE))
#' data$participant = 1
#' data$stimulus <- stimulus
#' data$discriminability <- discriminability
#' data <- data[data$response!=0, ] # drop not finished decision processes
#' data <- data[,-3] # drop conf measure (unobservable variable)
#' head(data)
#'
#' # 2. Define some parameter set in a data.frame
#' paramDf <- data.frame(a=2.5,v1=0.5, v2=1, t0=0.1,z=0.7,
#'                       sz=0,sv=0.2, st0=0,  tau=3, w=0.3,
#'                       theta1=0.8, svis=0.5, sigvis=0.8)
#'
#' # 3. Compute log likelihood for parameter and data
#' LogLikWEV(data, paramDf, model="dynWEV", condition="discriminability")
#' # adding the hypothetical interjudgment time to response times
#' # results in the same log likelihood as before when simult_conf=TRUE
#' data$rt <- data$rt + paramDf$tau
#' LogLikWEV(data, paramDf, model="dynWEV", condition="discriminability", simult_conf=TRUE)
#'
#' # the same function for "2DSD" model
#' paramDf <- data.frame(a=2.5,v1=0.5, v2=1, t0=0.1,z=0.7,
#'                       sz=0,sv=0.2, st0=0,  tau=3, theta1=0.8)
#' LogLikWEV(data, paramDf, model="2DSD", condition="discriminability", simult_conf=TRUE)
#' # this results in the same log likelihood as before
#'


#' @rdname LogLikWEV
#' @export
LogLikWEV <- function(data, paramDf, model="dynaViTE", simult_conf = FALSE, precision=1e-5, stop_on_error = TRUE, data_names = list(), ...) {
  #### Check data formatting ####
  data <- rename(data, ...)
  if ((model %in% c("dynWEV", "2DSD")) && !("lambda" %in% names(paramDf))) paramDf$lambda <- 0
  if (model=="dynWEV") model <- "dynaViTE"
  if (model=="2DSDT") model <- "2DSD"
  if (!("lambda" %in% names(paramDf))) paramDf$lambda <- 0


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
                          condition = as.numeric(factor(.data[[data_names$condition]],levels = sort(unique(data[[data_names$condition]])))))
  if (length(unique(data[[data_names$stimulus]]))==2) {
    data[[data_names$stimulus]] = if_else(data[[data_names$stimulus]]==sort(unique(data[[data_names$stimulus]]))[1],-1,1)
  }
  if (!is.numeric(data$rating)) {
    data <- data %>% mutate(rating = as.numeric(as.factor(.data[[data_names$rating]])))
  }
  ## Compute the row-wise likelihood of observations

  data <-data %>% mutate(vth1 = case_when(.data$response=="upper" ~ thetas_upper[.data$rating],
                                          .data$response=="lower" ~ thetas_lower[(.data$rating)]),
                         vth2 = case_when(.data$response=="upper" ~ thetas_upper[(.data$rating+1)],
                                          .data$response=="lower" ~ thetas_lower[(.data$rating+1)]),
                         M_drift = V[.data$condition]*.data$stimulus,
                         SV = SV[.data$condition],
                         S = S[.data$condition])


  probs <- with(data, switch(which(model== c("dynaViTE", "2DSD")),
                             dWEV(rt, vth1,vth2,
                                        response=response,
                                        tau=paramDf$tau, a=paramDf$a,
                                        v = M_drift,
                                        t0 = paramDf$t0, z = paramDf$z, sz = paramDf$sz, st0=paramDf$st0,
                                        sv = SV, w=paramDf$w, muvis=abs(M_drift), svis=paramDf$svis,
                                        sigvis=paramDf$sigvis, lambda=paramDf$lambda, s = S,
                                        simult_conf = simult_conf, z_absolute = FALSE,
                                        precision = precision, stop_on_error = stop_on_error,
                                        stop_on_zero=FALSE),
                             d2DSD(rt, response=response, vth1,vth2,
                                   tau=paramDf$tau, a=paramDf$a,
                                   v = M_drift,
                                   t0 = paramDf$t0, z = paramDf$z, sz = paramDf$sz, st0=paramDf$st0,
                                   sv = SV, lambda=paramDf$lambda, s=S,
                                   simult_conf = simult_conf, z_absolute = FALSE,
                                   precision = precision, stop_on_error = stop_on_error,
                                   stop_on_zero=FALSE)))



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
