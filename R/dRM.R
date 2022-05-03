#' Independent and partially anti-correlated Race Model for Decision Confidence
#'
#' Probability densities for response times, decisions and confidence judgments
#' in the independent Race Model (dIRM) or partially (anti-)correlated Race Model (dPCRM),
#' i.e. the probability of a given response (response: winning accumulator (1 or 2)) at a
#' given time (rts) and the confidence measure in interval between th1 and th2. The definition
#' of the confidence measure depends on the argument \code{time_scaled} (see Details).
#' The computations are based on Moreno-Bote (2010).
#'
#' The parameters for the models are \code{mu1} and \code{mu2} for the drift rates, \code{a}, \code{b}
#' for the starting points of the two accumulators and \code{s}/\code{sigma} for the incremental standard deviation
#' of the processes and \code{t0} and \code{st0} for the  minimum and range of uniformly distributed
#' non-decision times (including encoding and motor time).
#' For dIRM and dPCRM additionally the parameters \code{th1} and \code{th2} for
#' the lower and upper bound of the interval for the loosing accumulator.
#'
#' @param rt a numeric vector of RTs. For convenience also a \code{data.frame} with
#' columns \code{rt} and \code{response} (and \code{xj}) is possible.
#' @param response numeric vector of 1s and 2s, giving the accumulator that hits his boundary first.
#' @param mu1 numeric. Drift rate for the first accumulator
#' @param mu2 numeric. Drift rate for the second accumulator
#' @param a positive numeric. Upper threshold on the first accumulator.
#' @param b positive numeric. Upper threshold on the second accumulator.
#' @param th1 numeric. Lower bound of interval range for the confidence measure.
#' @param th2 numeric. Upper bound of interval range for the confidence measure.
#' @param wx numeric. Weight on losing accumulator for the computation of the confidence measure.
#'                    (Used only if time_scale==TRUE, Default = 1)
#' @param wrt numeric. Weight on reaction time for the computation of the confidence measure.
#'                    (Used only if time_scale==TRUE, Default wrt=0)
#' @param wint numeric. Weight on the interaction of losing accumulator and reaction time for
#'                      the computation of the confidence measure. (Used only if time_scale==TRUE, Default wint=0)
#' @param t0 numeric. Lower bound of non-decision time component in observable response times. Range: t0>=0. Default: t0=0.
#' @param st0 numeric. Range of a uniform distribution for non-decision time. Range: st0>=0. Default: st0=0.
#' @param s numeric. Diffusion constant of the two accumulators.  Usually fixed to 1 for most purposes as it scales other parameters (see Details). Range: s>0, Default: s=1.
#'
#' @param time_scaled logical. Whether the confidence measure should be time-dependent. See Details.
#' @param step_width numeric. Step size for the integration in t0 (motor time). Default: 1e-6.
#'
#' @return Returns the numerical value of the probability density in a numerical vector of the same length as
#' rt.
#'
#' @details \code{time_scaled} determines whether the confidence measure is computed in accordance to the
#' Balance of Evidence hypothesis, i.e. if response is 1 at time T, then \eqn{conf = b - X_2(T)}.
#' Otherwise, if \code{time_scaled=TRUE} (default), confidence is computed as linear combination of
#' Balance of Evidence and decision time, i.e.
#' \eqn{conf = wx (b-X_2 (T)) + wrt\frac{1}{\sqrt{T}} + wint\frac{b-X_2(T)}{\sqrt{T}}}.
#' Usually the weights (wx, wrt, wint) should sum to 1, as the confidence thresholds (th1 and th2) may be scaled
#' according to their sum. If this is not the case, they will be scaled accordingly internally! Usually, this
#' formula results in lower confidence when the reaction time is longer but the state of the second accumulator
#' is held constant. It is based on the optimal decision confidence in Moreno-Bote (2010).
#'
#' For convenience, the function allows that the first argument is a \code{data.frame} containing the
#' information of the first and second argument in the columns
#' (i.e., \code{rt} and \code{response}). Other columns (as well as passing
#' \code{response} separately as argument) will be ignored.
#'
#'
#' @note Similarly to the drift diffusion models (like \code{\link[rtdists:Diffusion]{ddiffusion}} and
#' \code{\link{dWEV}}), s is a scaling factor (scales: \code{mu1},\code{mu2}, \code{a},\code{b},
#' \code{th1},\code{th2},and \code{wrt}) and is usually fixed to 1.
#'
#' @references Moreno-Bote, R. (2010). Decision confidence and uncertainty in diffusion models with partially
#' correlated neuronal integrators. Neural Computation, 22(7), 1786–1811. https://doi.org/10.1162/neco.2010.12-08-930
#'
#'
#' @author R implementation and C code by Sebastian Hellmann.
#'
#' @useDynLib dynConfiR, .registration = TRUE
#'
#' @name dRM
#' @aliases dIRM dPCRM
#' @importFrom Rcpp evalCpp


#' @rdname dRM
#' @export
dIRM <- function (rt,response=1, mu1, mu2, a, b,
                  th1, th2, wx=1, wrt=0, wint=0,
                  t0=0, st0=0, s=1,
                  time_scaled = TRUE, step_width=NULL)
{
  # for convenience accept data.frame as first argument.
  if (is.data.frame(rt)) {
    response <- rt$response
    rt <- rt$rt
  }
  if (is.null(step_width)) {
    step_width = 0.089045 * exp(-1.037580*4)
  } else if (step_width>1) {
    step_width = 0.089045 * exp(-1.037580*step_width)
  }
  if (any(c(a<=0, b<=0))) {stop("Both thresholds (a  and b) must be positive")}
  if (any(s<=0)) {stop("s must be positive")}
  if (any(t0<0)) {stop("Non-decision time, t0, has to be non-negative")}
  if (any(st0<0)) {stop("Non-decision time range, st0, has to be non-negative")}


  nn <- length(rt)
  if (!time_scaled) {
    wrt <- 0
    wint <- 0
    wx <- 1
  }

  pars <- prepare_RaceModel_parameter(response = response,
                                      mu1, mu2,
                                      a, b,
                                      s, th1, th2,
                                      t0, st0,
                                      wx, wrt, wint,
                                      nn)
  densities <- vector("numeric",length=nn)
  for (i in seq_len(length(pars$parameter_indices))) {
    ok_rows <- pars$parameter_indices[[i]]

    densities[ok_rows] <- d_IRM (rt[ok_rows]-pars$params[ok_rows[1],12],
                                 pars$params[ok_rows[1],1:11],
                                 pars$params[ok_rows[1],13],
                                 step_width)
  }
  abs(densities)
}

#' @rdname dRM
#' @export
dPCRM <- function (rt,response=1, mu1, mu2, a, b,
                   th1, th2, wx=1, wrt=0, wint=0,
                   t0=0, st0=0, s=1,
                   time_scaled = TRUE, step_width=NULL)
{
  # for convenience accept data.frame as first argument.
  if (is.data.frame(rt)) {
    response <- rt$response
    rt <- rt$rt
  }
  if (is.null(step_width)) {
    step_width = 0.089045 * exp(-1.037580*4)
  } else if (step_width>1) {
    step_width = 0.089045 * exp(-1.037580*step_width)
  }
  if (any(c(a<=0, b<=0))) {stop("Both thresholds (a  and b) must be positive")}
  if (any(s<=0)) {stop("s must be positive")}
  if (any(t0<0)) {stop("Non-decision time, t0, has to be non-negative")}
  if (any(st0<0)) {stop("Non-decision time range, st0, has to be non-negative")}

  nn <- length(rt)
  if (!time_scaled) {
    wrt <- 0
    wint <- 0
    wx <- 1
  }

  pars <- prepare_RaceModel_parameter(response = response,
                                      mu1, mu2,
                                      a, b,
                                      s, th1, th2,
                                      t0, st0,
                                      wx, wrt, wint,
                                      nn)
  densities <- vector("numeric",length=nn)
  for (i in seq_len(length(pars$parameter_indices))) {
    ok_rows <- pars$parameter_indices[[i]]

    densities[ok_rows] <- d_PCRM (rt[ok_rows]-pars$params[ok_rows[1],12],
                                  pars$params[ok_rows[1],1:11],
                                  pars$params[ok_rows[1],13],
                                  step_width)
  }
  abs(densities)
}

