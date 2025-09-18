#' Correlated Multiple-threshold Log-normal Race Model for Decision Confidence
#'
#' Probability densities and random number generators for response times,
#' decisions and confidence judgments in the multi-threshold correlated
#' log-normal race model (`MTLNR`; Reynolds et al., 2020),
#' i.e. the probability of a given response (response: winning accumulator
#' (1 or 2)) at a given time (rt) and the confidence measure in the interval
#' between `th1` and `th2`. The confidence measure is defined as the log-ratio between
#' the time the loosing accumulator would hit its boundary and the decision time
#' (which means it is in the interval (0, Inf).
#' The parameters for the model are:
#' \code{mu_v1}/\code{mu_v2} and \code{s_v1}/\code{s_v2} the mean and standard deviation
#' parameters for the log-normally distributed accumulation rates of the two accumulators,
#' \code{mu_d1}/\code{mu_d2} and \code{s_d1}/\code{s_d2} the mean and standard deviation
#' parameters for the log-normally distributed boundary distances of the two accumulators,
#' \code{rho_v} and \code{rho_d} the correlation coefficients for the accumulation rates and
#' the boundary distance, giving the correlation of these parameters between the two accumulators,
#' \code{t0} and \code{st0} for the minimum and range of uniformly distributed
#' non-decision times (including encoding and motor time).
#'
#' For the computation of confidence judgments, the parameters for the confidence
#' thresholds \code{th1} and
#' \code{th2} for the lower and upper bound of the interval for the confidence
#' measure, or `thresholds` for a vector of thresholds (see Details).
#'
#' @param rt a numeric vector of RTs. For convenience also a \code{data.frame} with
#' columns \code{rt} and \code{response} is possible.
#' @param response numeric vector with values in `c(1, 2)`, giving the accumulator that hit its
#' boundary first.
#' @param rating numeric vector with integer values from 1 to the number of confidence levels
#' @param th1 numeric. Lower bound of interval range for the confidence measure.
#' @param th2 numeric. Upper bound of interval range for the confidence measure.
#' @param thresholds numeric vector of length 2*(number of confidence levels).
#' Confidence thresholds, which are used to compare the confidence variable against
#' for producing discrete confidence judgments. The first half of the entries should
#' be increasing and represent the confidence thresholds for a `response` of 1; the second
#' half of the entries should also be increasing and are the thresholds for `response`=2.
#' @param mu_v1 numeric. Mean parameter of log-normally distributed accumulation rate of first accumulator
#' @param mu_v2 numeric. Mean parameter of log-normally distributed accumulation rate of second accumulator
#' @param s_v1  numeric. Standard deviation of log-normally distributed accumulation rate of first accumulator
#' @param s_v2  numeric. Standard deviation of log-normally distributed accumulation rate of second accumulator
#' @param rho_v numeric. Correlation parameter for accumulation rates
#' @param mu_d1 numeric. Mean parameter of log-normally distributed boundary distance of first accumulator
#' @param mu_d2 numeric. Mean parameter of log-normally distributed boundary distance of second accumulator
#' @param s_d1 numeric.  Standard deviation of log-normally distributed boundary distance of first accumulator
#' @param s_d2 numeric.  Standard deviation of log-normally distributed boundary distance of second accumulator
#' @param rho_d numeric. Correlation parameter for boundary distances
#' @param t0 numeric. Lower bound of non-decision time component in observable response times.
#' Range: `t0>=0`. Default: 0.
#' @param st0 numeric. Range of a uniform distribution for non-decision time. Range: `st0>=0`.
#' Default: 0.
#' @param precision numerical scalar value. Precision of calculation. Determines the
#' step size of integration w.r.t. `t0`.
#' Represents the number of decimals precisely computed on average. Default is 6.
#' @param step_width numeric. Alternative way to define the precision of integration
#' w.r.t. `t0` by directly providing the step size for the integration.
#'
#' @param n integer. The number of samples generated.
#'
#' @return `dMTLNR` and `dMTLNR_multiple_ratings` return the numerical value of the probability density in a numerical vector of the same
#' length as `rt`.
#'
#' `rMTLNR` returns a `data.frame` with five columns and `n` rows. Column names are `rt` (response
#' time), `response` (1 or 2, indicating which accumulator hit its boundary first),
#' `Tdec` (the actual decision time (without non-decision time),
#' `conf` (the log of the ratio of boundary hitting times (the confidence variable),
#' and `rating` (the discrete confidence judgment).
#'
#' The race parameters (as well as \code{response} (and \code{rating}), \code{th1},
#' and \code{th2}) are recycled to the length of the result (either `rt` or `n`).
#' In other words, the functions are completely vectorized for all parameters
#' and even the response.
#'
#' @details
#' The model assumes that each of the two accumulators has a log-normally distributed
#' boundary distance \eqn{D}{D} with mean parameter \eqn{\mu_D}{mu_D} and standard
#' deviation parameter \eqn{s_D^2}{s_D^2} and a log-normally distributed accumulation rate
#' \eqn{V}{V} with respective parameters. The accumulation takes the form of a linear
#' ballistic accumulation without any noise, such that the boundary crossing times
#' \eqn{T=D/V}{T=D/V} are log-normally distributed with mean parameter
#' \eqn{\mu_D - \mu_V}{mu_D-mu_V} and variance parameter
#' \eqn{s_D^2 + s_V^2}{s_D^2 + s_V^2}.
#' In addition, the boundary distances for the two accumulators are correlated with
#' the correlation determined by \eqn{\rho_D}{rho_D}. Similarly, the accumulation
#' rates share a correlation with parameter \eqn{\rho_V}{rho_V}.
#' Confidence is determined by the log-ratio of the loosing over the winning boundary
#' crossing time, i.e., if the first accumulator hit its boundary first, confidence
#' is determined by
#' \deqn{conf = log(T_2 / T_1).}{conf = log(T_2 / T_1).}
#' This confidence measure is then compared to the set of `thresholds` to produce
#' discrete confidence judgments.
#' This is equivalent to a confident computation based on the Balance of Evidence
#' at decision time, although symmetry conditions for the thresholds may differ
#' depending on the interpretation (see Reynolds et al., 2020 for more detail).
#'
#' For convenience, the likelihood function allows that the first argument is a \code{data.frame} containing the
#' information of the first and second argument in the columns
#' (i.e., \code{rt} and \code{response} (and \code{rating} if relevant)).
#' Other columns (as well as passing \code{response} separately as argument)
#' will be ignored.
#'
#' \emph{Difference between `dMTLNR` and `dMTLNR_multiple_ratings`}
#' The function `dMTLNR` allows to compute the probability of a rt and response with
#' the confidence variable being within an interval given by two thresholds, `th1`
#' and `th2`, similar to the definitions of the other density functions
#' (like \code{\link{ddynaViTE}}).
#' The function `dMTLNR_multiple_ratings` takes a vector for the discrete confidence
#' judgments, `rating`, and a vector `thresholds`, such that the confidence
#' interval can vary from observation to observation. The correct interval limits
#' are chosen by the function depending on the entry in `rating`.
#'
#' @note The model is highly over-parametrized because the mean parameters for
#' the boundary distances and accumulation rates trade off. Similarly, the
#' variance parameters and correlation parameters trade off. For this reason,
#' one may only use the first set of parameters for the accumulation rates (`mu_v1`,...).
#'
#' @references  Reynolds, A., Kvam, P. D., Osth, A. F., & Heathcote, A. (2020). Correlated racing evidence accumulator models. \emph{Journal of Mathematical Psychology, 96}, 102331. doi: doi: 10.1016/j.jmp.2020.102331
#'
#' @author Sebastian Hellmann
#'
#' @useDynLib dynConfiR, .registration = TRUE
#'
#' @name MTLNR
#' @aliases lognormalrace
#' @importFrom Rcpp evalCpp
#'
#' @examples
#' # Plot rt distribution ignoring confidence
#' curve(dMTLNR(x, 1, th1=-Inf, th2=Inf, mu_v1=0.5, mu_v2=-0.5,
#'             mu_d1=1, mu_d2=1, t0=0.1), xlim=c(0,2.5), ylim=c(0,2))
#' curve(dMTLNR(x, 2, th1=-Inf, th2=Inf, mu_v1=0.5, mu_v2=-0.5,
#'             mu_d1=1, mu_d2=1, t0=0.1), col="red", add=TRUE)
#' # t0 indicates minimal response time possible
#' abline(v=0.1)
#'
#' # Generate a random sample
#' df1 <- rMTLNR(5000, mu_v1=0.2, mu_v2=-0.2, mu_d1=1, mu_d2=1, t0=0.1)
#' head(df1)
#'
#' # Compute density with rt and response as separate arguments
#' dMTLNR(seq(0, 2, by =0.4), response=2, th1=0.5, th2=2,
#'       mu_v1=0.2, mu_v2=-0.2, mu_d1=1, mu_d2=1, t0=0.1)
#'
#' # Compute density with rt and response in data.frame argument
#' df1 <- subset(df1, response !=0) # drop trials where no accumulation hit its boundary
#' dMTLNR(df1[1:5,], th1=0, th2=Inf, mu_v1=0.2, mu_v2=-0.2,
#'       mu_d1=1, mu_d2=1, t0=0.1)
#'
#' # Example with correlation parameters
#' dMTLNR(df1[1:5,], th1=0, th2=Inf, mu_v1=0.2, mu_v2=-0.2,
#'       mu_d1=1, mu_d2=1, rho_v=0.3, rho_d=0.2, t0=0.1)
#'
#' # Example with multiple confidence ratings using dMTLNR_multiple_ratings
#' thresholds <- c(0.5, 1.5, 2.5,  # for response=1 (increasing)
#'                0.3, 1.2, 2.0)  # for response=2 (increasing)
#'
#' # Create some sample data with ratings
#' sample_data <- data.frame(
#'  rt = c(0.8, 1.2, 0.9, 1.5, 1.1),
#'  response = c(1, 2, 1, 2, 1),
#'  rating = c(1, 2, 3, 1, 2)
#' )
#'
#' dMTLNR_multiple_ratings(sample_data, thresholds=thresholds,
#'                        mu_v1=0.2, mu_v2=-0.2, mu_d1=1, mu_d2=1, t0=0.1)
#'
#' # Compare RT and confidence distributions for different parameter settings
#' df_low_var <- rMTLNR(2000, thresholds = thresholds,
#'                     mu_v1=0.3, mu_v2=-0.3, mu_d1=1, mu_d2=1,
#'                     s_v1=0.5, s_v2=0.5, s_d1=0.5, s_d2=0.5, t0=0.1)
#' df_high_var <- rMTLNR(2000, thresholds=thresholds,
#'                      mu_v1=0.3, mu_v2=-0.3, mu_d1=1, mu_d2=1,
#'                      s_v1=1.5, s_v2=1.5, s_d1=1.5, s_d2=1.5, t0=0.1)
#'
#' two_samples <- rbind(cbind(df_low_var, variance="low"),
#'                     cbind(df_high_var, variance="high"))
#' two_samples <- two_samples[two_samples$response != 0, ]
#'
#' # Compare RT distributions
#' boxplot(log(rt) ~ variance + response, data = two_samples)
#'
#' # Compare confidence distributions
#' boxplot(conf ~ variance + response, data = two_samples)
#' boxplot(rating ~ variance + response, data = two_samples)
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'  require(ggplot2)
#'  ggplot(two_samples, aes(x = rt, y = conf)) +
#'   stat_density_2d(aes(fill = after_stat(density)), geom = "raster", contour = FALSE, na.rm=TRUE) +
#'   facet_grid(cols = vars(variance), rows = vars(response),
#'              labeller = "label_both") +
#'   xlim(c(0.2, 2.0)) + ylim(c(-2, 4))
#' }
#'

#' @rdname MTLNR
#' @export
dMTLNR <- function (rt,response=1,
                    th1, th2,
                    mu_v1, mu_v2,
                    s_v1=1, s_v2=1, rho_v=0,
                    mu_d1=0, mu_d2=0, s_d1=1, s_d2=1,
                    rho_d=0,
                    t0=0, st0=0,
                    precision=6, step_width=NULL)
{
  # for convenience accept data.frame as first argument.
  if (is.data.frame(rt)) {
    response <- rt$response
    rt <- rt$rt
  }
  if (length(response)!=length(rt)) {
    if (length(response)!=1) warning("Length of response is neither 1 nor length of rt, and will be lep_len'ed")
    response = rep_len(response, length(rt))
  }

  if (is.null(step_width)) {
    step_width = 0.4578627727708822 * exp(-0.8466054539598147*precision)
    #step_width = 0.089045 * exp(-1.037580*precision)
    # step_width = 0.089045 * exp(-1.037580*4)
  } else if (step_width>1) {
    step_width = 0.4578627727708822 * exp(-0.8466054539598147*step_width)
    #step_width = 0.089045 * exp(-1.037580*step_width)
  }

  # print(step_width)
  if (any(c(s_v1<0, s_v2<0))) {stop("Variances of all accumulation rates must be positive")}
  if (any(c(s_d1<0, s_d2<0))) {stop("Variances of all boundary distances must be positive")}
  if (all(c(s_d1, s_d2, s_v1, s_v2)==0)) {stop("Not all variance components can be 0!")}
  if (any(c(rho_d < (-1),rho_d > 1) )) {stop("Correlation of boundary distances must be in the interval -1 and 1")}
  if (any(c(rho_v < (-1),rho_v > 1) )) {stop("Correlation of accumulation rates must be in the interval -1 and 1")}
  if (any(t0<0)) {stop("Non-decision time, t0, has to be non-negative")}
  if (any(st0<0)) {stop("Non-decision time range, st0, has to be non-negative")}

  nn <- length(rt)
  pars <- prepare_MTLNR_parameter(response = response,
                                  th1, th2,
                                  mu_v1, mu_v2,
                                  s_v1, s_v2, rho_v,
                                  mu_d1, mu_d2, s_d1, s_d2,
                                  rho_d,
                                  t0, st0,
                                  nn)
  densities <- vector("numeric",length=nn)
  for (i in seq_len(length(pars$parameter_indices))) {
    ok_rows <- pars$parameter_indices[[i]]

    densities[ok_rows] <- d_MT_LNR (rt[ok_rows],rating=rep(1, length(ok_rows)),
                                    params=pars$params[ok_rows[1],1:12],
                                    thresholds = pars$params[ok_rows[1],14:15],
                                 win=pars$params[ok_rows[1],13],
                                 step_width)
  }
  abs(densities)
}

#' @rdname MTLNR
#' @export
dMTLNR_multiple_ratings <-
  function (rt,response=1, rating=1, thresholds,
                    mu_v1, mu_v2,
                    s_v1=1, s_v2=1, rho_v=0,
                    mu_d1=0, mu_d2=0, s_d1=1, s_d2=1,
                    rho_d=0,
                    t0=0, st0=0,
                    precision=6, step_width=NULL)
{
  # for convenience accept data.frame as first argument.
  if (is.data.frame(rt)) {
    response <- rt$response
    rating <- rt$rating
    rt <- rt$rt
  }

  if (!all(rating==as.integer(rating))) {
    warning("Entries in rating will be converted to integer values!")
  }
  rating <- as.integer(rating)
  if (length(rating)!=length(rt)) {
    if (length(rating)!=1) warning("Length of rating is neither 1 nor length of rt, and will be lep_len'ed")
    rating = rep_len(rating, length(rt))
  }
  if (length(response)!=length(rt)) {
    if (length(response)!=1) warning("Length of response is neither 1 nor length of rt, and will be lep_len'ed")
    response = rep_len(response, length(rt))
  }

  if (is.null(step_width)) {
    step_width = 0.4578627727708822 * exp(-0.8466054539598147*precision)
    #step_width = 0.089045 * exp(-1.037580*precision)
    # step_width = 0.089045 * exp(-1.037580*4)
  } else if (step_width>1) {
    step_width = 0.4578627727708822 * exp(-0.8466054539598147*step_width)
    #step_width = 0.089045 * exp(-1.037580*step_width)
  }

  # print(step_width)
  if (any(c(s_v1<0, s_v2<0))) {stop("Variances of all accumulation rates must be positive")}
  if (any(c(s_d1<0, s_d2<0))) {stop("Variances of all boundary distances must be positive")}
  if (all(c(s_d1, s_d2, s_v1, s_v2)==0)) {stop("Not all variance components can be 0!")}
  if (any(c(rho_d < (-1),rho_d > 1) )) {stop("Correlation of boundary distances must be in the interval -1 and 1")}
  if (any(c(rho_v < (-1),rho_v > 1) )) {stop("Correlation of accumulation rates must be in the interval -1 and 1")}
  if (any(t0<0)) {stop("Non-decision time, t0, has to be non-negative")}
  if (any(st0<0)) {stop("Non-decision time range, st0, has to be non-negative")}
  if (any(thresholds < 0)) { stop("thresholds need to be greater than or equal 0")}
  if (any(diff(thresholds[1:(length(thresholds)/2)]) < 0)) { stop("thresholds need to be in ascending order (within each response); first half of thresholds (for 1 responses) are not")}
  if (any(diff(thresholds[(length(thresholds)/2+1):length(thresholds)]) < 0)) { stop("thresholds need to be in ascending order (within each response); second half of thresholds (for 2 responses) are not")}

  if (!all(c(length(mu_v1) ,
             length(mu_v2) ,
             length(s_v1 ) ,
             length(s_v2 ) ,
             length(rho_v) ,
             length(mu_d1) ,
             length(mu_d2) ,
             length(s_d1 ) ,
             length(s_d2 ) ,
             length(rho_d) ,
             length(t0   ) ,
             length(st0))==1)) stop("For this function to work, all parameters (except for thresholds) have to have length 1")

  if (is.null(thresholds)) {
    warning("thresholds was NULL, so we continue assuming no thresholds and all confidence judgments are always 1.")
    thresholds_1 <- thresholds_2 <- c(0, .Machine$double.xmax)
    if (!all(rating==1)) stop("Confidence is not always 1, but no confidence thresholds provided!")
  } else {
    thresholds_1 <- c(0, thresholds[1:(length(thresholds)/2)], .Machine$double.xmax)
    thresholds_2 <- c(0, thresholds[(length(thresholds)/2+1):length(thresholds)], .Machine$double.xmax)
  }

  nn <- length(rt)
  densities <- vector("numeric",length=nn)

  ind_win1  <- response==1
    densities[ind_win1] <- d_MT_LNR(rt[ind_win1], rating[ind_win1],
                                  params=c(mu_d1, mu_d2, mu_v1, mu_v2,
                                           s_d1 , s_d2 , s_v1 , s_v2 ,
                                           rho_d, rho_v, t0   , st0   ),
                                  thresholds_1,
                                 win=1,
                                 step_width)
  densities[!ind_win1] <- d_MT_LNR(rt[!ind_win1], rating[!ind_win1],
                                   params=c(mu_d1, mu_d2, mu_v1, mu_v2,
                                            s_d1 , s_d2 , s_v1 , s_v2 ,
                                            rho_d, rho_v, t0   , st0   ),
                                   thresholds_2,
                                  win=2,
                                  step_width)

  return(abs(densities))
}





#' @rdname MTLNR
#' @export
rMTLNR <- function (n, thresholds=NULL,
                    mu_v1, mu_v2,
                    s_v1=1, s_v2=1, rho_v=0,
                    mu_d1=0, mu_d2=0, s_d1=1, s_d2=1,
                    rho_d=0,
                    t0=0, st0=0)
{
  if (any(missing(mu_v1), missing(mu_v2))) stop("mu_v1 and mu_v2 must be supplied")
  if (any(c(s_v1<0, s_v2<0))) {stop("Variances of all accumulation rates must be positive")}
  if (any(c(s_d1<0, s_d2<0))) {stop("Variances of all boundary distances must be positive")}
  if (all(c(s_d1, s_d2, s_v1, s_v2)==0)) {stop("Not all variance components can be 0!")}
  if (any(c(rho_d < (-1),rho_d > 1) )) {stop("Correlation of boundary distances must be in the interval -1 and 1")}
  if (any(c(rho_v < (-1),rho_v > 1) )) {stop("Correlation of accumulation rates must be in the interval -1 and 1")}
  if (any(t0<0)) {stop("Non-decision time, t0, has to be non-negative")}
  if (any(st0<0)) {stop("Non-decision time range, st0, has to be non-negative")}

  if (missing(thresholds)|| is.null(thresholds)) thresholds <- 1e+128
  pars <- prepare_MTLNR_parameter(response = 1,
                                       1, 99,
                                       mu_v1, mu_v2,
                                       s_v1, s_v2, rho_v,
                                       mu_d1, mu_d2, s_d1, s_d2,
                                       rho_d,
                                       t0, st0,
                                      n)
  res <- matrix(NA, nrow=n, ncol=5)
  for (i in seq_len(length(pars$parameter_indices))) {
    ok_rows <- pars$parameter_indices[[i]]
    current_n <- length(ok_rows)
    out <- r_MT_LNR(current_n, pars$params[ok_rows[1],c(1:12)],
                    thresholds)
    res[ok_rows,1:5] <- out

  }

  res <- as.data.frame(res)
  names(res) <- c("rt", "response", "Tdec", "conf", "rating")


  return(res)
}





prepare_MTLNR_parameter <- function(response,th1, th2,
                                         mu_v1, mu_v2,
                                         s_v1, s_v2, rho_v,
                                         mu_d1, mu_d2, s_d1, s_d2,
                                         rho_d,
                                         t0, st0, nn) {

  if ((length(mu_v1) == 1) &
      (length(mu_v2) == 1) &
      (length(s_v1 ) == 1) &
      (length(s_v2 ) == 1) &
      (length(rho_v) == 1) &
      (length(mu_d1) == 1) &
      (length(mu_d2) == 1) &
      (length(s_d1 ) == 1) &
      (length(s_d2 ) == 1) &
      (length(rho_d) == 1) &
      (length(t0   ) == 1) &
      (length(st0  ) == 1) &
      (length(th1) ==1)    &
      (length(th2) ==1) ) {
    skip_checks <- TRUE
  } else {
    skip_checks <- FALSE
  }

  response <- as.numeric(response)
  if (any(!(response %in% 1:2)))
    stop("response needs to be  %in% 1:2!")
  numeric_bounds <- as.integer(response)

  numeric_bounds <- rep(numeric_bounds, length.out = nn)
  if (!skip_checks) {
    # all parameters brought to length of rt
    mu_v1 <- rep(mu_v1, length.out=nn)
    mu_v2 <- rep(mu_v2, length.out=nn)
    s_v1  <- rep(s_v1 , length.out=nn)
    s_v2  <- rep(s_v2 , length.out=nn)
    rho_v <- rep(rho_v, length.out=nn)
    mu_d1 <- rep(mu_d1, length.out=nn)
    mu_d2 <- rep(mu_d2, length.out=nn)
    s_d1  <- rep(s_d1 , length.out=nn)
    s_d2  <- rep(s_d2 , length.out=nn)
    rho_d <- rep(rho_d, length.out=nn)
    t0    <- rep(t0   , length.out=nn)
    st0   <- rep(st0  , length.out=nn)
    th1   <- rep(th1  , length.out = nn)
    th2   <- rep(th2  , length.out = nn)
  }
  th1[th1==-Inf] <- 0
  th2[th2==Inf] <- .Machine$double.xmax
  # Build parameter matrix
  params <- cbind (mu_d1, mu_d2,
    mu_v1, mu_v2,
    s_d1 , s_d2 ,
    s_v1 , s_v2 ,
    rho_d, rho_v,
    t0   , st0,
    numeric_bounds, th1, th2)

  # Check for illegal parameter values
  if(ncol(params)<15) stop("Not enough parameters supplied: probable attempt to pass NULL values?")
  if(!is.numeric(params)) stop("Parameters need to be numeric.")
  if (any(is.na(params)) || !all(is.finite(params))) {
    stop("Parameters need to be numeric and finite.")
  }


  if (!skip_checks) {
    parameter_char <- apply(params, 1, paste0, collapse = "\t")
    parameter_factor <- factor(parameter_char, levels = unique(parameter_char))
    parameter_indices <- split(seq_len(nn), f = parameter_factor)
  } else {
    if (all(numeric_bounds == 2L) | all(numeric_bounds == 1L)) {
      parameter_indices <- list(
        seq_len(nn)
      )
    } else {
      parameter_indices <- list(
        seq_len(nn)[numeric_bounds == 2L],
        seq_len(nn)[numeric_bounds == 1L]
      )
    }
  }
  list(
    params = params
    , parameter_indices = parameter_indices
  )
}





