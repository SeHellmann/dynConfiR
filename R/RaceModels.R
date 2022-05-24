#' Independent and partially anti-correlated Race Model for Decision Confidence
#'
#' Probability densities and random number generators for response times,
#' decisions and confidence judgments in the independent Race Model
#' (`dIRM`/`rIRM`) or partially (anti-)correlated Race Model (`dPCRM`/`rPCRM`),
#' i.e. the probability of a given response (response: winning accumulator
#' (1 or 2)) at a given time (rt) and the confidence measure in the interval
#' between th1 and th2 (Hellmann et al., preprint). The definition of the
#' confidence measure depends on the argument \code{time_scaled} (see Details).
#' The computations are based on Moreno-Bote (2010).
#' The parameters for the models are \code{mu1} and \code{mu2} for the drift
#' rates, \code{a}, \code{b} for the upper thresholds of the two accumulators
#' and \code{s} for the incremental standard deviation of the processes and
#' \code{t0} and \code{st0} for the  minimum and range of uniformly distributed
#' non-decision times (including encoding and motor time).
#' For the computation of confidence judgments, the parameters \code{th1} and
#' \code{th2} for the lower and upper bound of the interval for confidence
#' measure and if `time_scaled` is `TRUE` the weight parameters `wx`, `wrt`,
#' `wint` for the computation of the confidence measure are required (see Details).
#'
#' @param rt a numeric vector of RTs. For convenience also a \code{data.frame} with
#' columns \code{rt} and \code{response} is possible.
#' @param response numeric vector with values in `c(1, 2)`, giving the accumulator that hit its
#' boundary first.
#' @param mu1 numeric. Drift rate for the first accumulator
#' @param mu2 numeric. Drift rate for the second accumulator
#' @param a positive numeric. Distance from starting point to boundary of the first accumulator.
#' @param b positive numeric. Distance from starting point to boundary of the second accumulator.
#' @param th1 numeric. Lower bound of interval range for the confidence measure.
#' @param th2 numeric. Upper bound of interval range for the confidence measure.
#' @param wx numeric. Weight on losing accumulator for the computation of the confidence measure.
#'                    (Used only if time_scale==TRUE, Default = 1)
#' @param wrt numeric. Weight on reaction time for the computation of the confidence measure.
#'                    (Used only if time_scale==TRUE, Default wrt=0)
#' @param wint numeric. Weight on the interaction of losing accumulator and reaction time for
#'                      the computation of the confidence measure. (Used only if time_scale==TRUE,
#'                      Default wint=0)
#' @param t0 numeric. Lower bound of non-decision time component in observable response times.
#' Range: t0>=0. Default: t0=0.
#' @param st0 numeric. Range of a uniform distribution for non-decision time. Range: st0>=0.
#' Default: st0=0.
#' @param s numeric. Diffusion constant of the two accumulators.  Usually fixed to 1 for most
#' purposes as it scales other parameters (see Details). Range: s>0, Default: s=1.
#'
#' @param time_scaled logical. Whether the confidence measure should be time-dependent. See Details.
#' @param step_width numeric. Step size for the integration in t0 (motor time). Default: 1e-6.
#'
#' @param n integer. The number of samples generated.
#' @param delta numeric. Discretization step size for simulations in the stochastic process
#' @param maxrt numeric. Maximum decision time returned. If the simulation of the stochastic
#' process exceeds a decision time of `maxrt`, the `response` will be set to 0 and the `maxrt`
#' will be returned as `rt`.
#'
#' @return `dIRM` and `dPCRM` return the numerical value of the probability density in a numerical vector of the same
#' length as `rt`.
#'
#' `rIRM` and `dPCRM` return a `data.frame` with four columns and `n` rows. Column names are `rt` (response
#' time), `response` (1 or 2, indicating which accumulator hit its boundary first),
#' `xl` (the final state of the loosing accumulator), and `conf` (the
#' value of the confidence measure; not discretized!).
#'
#' The race parameters (as well as \code{response}, \code{th1},
#' and \code{th2}) are recycled to the length of the result (either `rt` or `n`).
#' In other words, the functions are completely vectorized for all parameters
#' and even the response.
#'
#' @details The parameters are formulated, s.t. both accumulators start at 0 and trigger a decision at their
#' positive boundary `a` and `b` respectively. That means, both parameters have to be positive.
#' Internally the computations adapt the parametrization of Moreno-Bote (2010).
#'
#' \code{time_scaled} determines whether the confidence measure is computed in accordance to the
#' Balance of Evidence hypothesis (if `time_scaled=FALSE`), i.e. if `response` is 1 at time T and
#' \eqn{X_2}{X2} is the second accumulator, then
#' \deqn{conf = b - X_2(T)}{conf = b - X2(T)}.
#' Otherwise, if \code{time_scaled=TRUE} (default), confidence is computed as linear combination of
#' Balance of Evidence, decision time, and an interaction term, i.e.
#' \deqn{conf = wx (b-X_2 (T)) + wrt\frac{1}{\sqrt{T}} + wint\frac{b-X_2(T)}{\sqrt{T}}.}{conf = wx (b-X2(T)) + wrt/\sqrt{T} + wint (b-X2(T))\sqrt{T}.}
#' Usually the weights (`wx`, `wrt`, `wint`) should sum to 1, as the confidence thresholds
#' (`th1` and `th2`) may be scaled according to their sum. If this is not the case, they will be scaled
#' accordingly internally! Usually, this formula results in lower confidence when the reaction time is
#' longer but the state of the second accumulator is held constant. It is based on the optimal decision
#' confidence in Moreno-Bote (2010).
#'
#'
#' For convenience, the likelihood function allows that the first argument is a \code{data.frame} containing the
#' information of the first and second argument in the columns
#' (i.e., \code{rt} and \code{response}). Other columns (as well as passing
#' \code{response} separately as argument) will be ignored.
#'
#' The simulations are done by simulating normal variables in discretized steps until
#' one process reaches the boundary. If no boundary is met within the maximum time, response is
#' set to 0.
#'
#' @note Similarly to the drift diffusion models (like \code{ddiffusion} and
#' \code{\link{dWEV}}), s is a scaling factor (scales: \code{mu1},\code{mu2}, \code{a},\code{b},
#' \code{th1},\code{th2},and \code{wrt}) and is usually fixed to 1.
#'
#' @references Hellmann, S., Zehetleitner, M., & Rausch, M. (preprint). Simultaneous modeling of choice,
#' confidence and response time in visual perception. https://osf.io/9jfqr/
#'
#' Moreno-Bote, R. (2010). Decision confidence and uncertainty in diffusion models with partially
#' correlated neuronal integrators. Neural Computation, 22(7), 1786–1811. https://doi.org/10.1162/neco.2010.12-08-930
#'
#'
#' @author Sebastian Hellmann
#'
#' @useDynLib dynConfiR, .registration = TRUE
#'
#' @name RaceModels
#' @aliases dIRM dPCRM rRM dRM rIRM rPCRM racemodels
#' @importFrom Rcpp evalCpp
#'
#' @examples
#' # Plot rt distribution ignoring confidence
#' curve(dPCRM(x, 1, mu1=0.5, mu2=-0.5, a=1, b=1, th1=-Inf, th2=Inf, t0=0.1), xlim=c(0,2.5))
#' curve(dPCRM(x, 2, mu1=0.5, mu2=-0.5, a=1, b=1, th1=-Inf, th2=Inf, t0=0.1), col="red", add=TRUE)
#' curve(dIRM(x, 1, mu1=0.5, mu2=-0.5, a=1, b=1, th1=-Inf, th2=Inf, t0=0.1), lty=2,add=TRUE)
#' curve(dIRM(x, 2, mu1=0.5, mu2=-0.5, a=1, b=1, th1=-Inf, th2=Inf, t0=0.1),
#'       col="red", lty=2, add=TRUE)
#' # t0 indicates minimal response time possible
#' abline(v=0.1)
#' ## Following example may be equivalently used for the IRM model functions.
#' # Generate a random sample
#' df1 <- rPCRM(5000,  mu1=0.2, mu2=-0.2, a=1, b=1, t0=0.1,
#'             wx = 1) # Balance of Evidence
#' # Same RT and response distribution but different confidence distribution
#' df2 <- rPCRM(5000,  mu1=0.2, mu2=-0.2, a=1, b=1, t0=0.1,
#'              wint = 0.2, wrt=0.8)
#' head(df1)
#'
#' # s scales other decision relevant parameters
#' dPCRM(df1[1:5,], 2, mu1=0.2, mu2=-0.2, a=1, b=1, th1=0, th2=Inf, t0=0.1)
#' s <- 2
#' dPCRM(df1[1:5,], 2, mu1=0.2*s, mu2=-0.2*s, a=1*s, b=1*s, th1=0, th2=Inf, t0=0.1, s=s)
#' # s also scales confidence parameters
#' dPCRM(df1[1:5,], 2, mu1=0.2, mu2=-0.2, a=1, b=1,
#'      th1=0.5, th2=2, wx = 0.3, wint=0.4, wrt=0.1, t0=0.1)
#' s <- 2
#' dPCRM(df1[1:5,], 2, mu1=0.2*s, mu2=-0.2*s, a=1*s, b=1*s,
#'      th1=0.5*s, th2=2*s, wx = 0.3, wint=0.4, wrt=0.1*s, t0=0.1, s=s)
#'
#' two_samples <- rbind(cbind(df1, ws="BoE"),
#'                    cbind(df2, ws="RT"))
#' # drop not finished decision processes
#' two_samples <- two_samples[two_samples$response!=0,]
#' # no difference in RT distributions
#' boxplot(rt~ws+response, data=two_samples)
#' # but different confidence distributions
#' boxplot(conf~ws+response, data=two_samples)
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#' require(ggplot2)
#' ggplot(two_samples, aes(x=rt, y=conf))+
#'   stat_density_2d(aes(fill = ..density..), geom = "raster", contour = FALSE, h=c(0.3, 0.7)) +
#'   xlim(c(0.2, 1.3))+ ylim(c(0, 2.5))+
#'   facet_grid(cols=vars(ws), rows=vars(response), labeller = "label_both")
#' }
#' # Restricting to specific confidence region
#' df1 <- df1[df1$conf >0 & df1$conf <1,]
#' dPCRM(df1[1:5,], th1=0, th2=1,mu1=0.2, mu2=-0.2, a=1, b=1, t0=0.1,wx = 1 )
#'


#' @rdname RaceModels
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

#' @rdname RaceModels
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


#' @rdname RaceModels
#' @export
rIRM <- function (n, mu1, mu2, a, b,
                 wx=1, wrt=0, wint=0,
                 t0=0, st0=0, s=1,
                 time_scaled = TRUE, step_width=NULL,
                 delta=0.01, maxrt=15)
{
  if (any(missing(mu1), missing(mu2),
          missing(a), missing(b))) stop("mu1, mu2, a, and b must be supplied")
  if (any(c(a<=0, b<=0))) {stop("Both thresholds (a  and b) must be positive")}
  if (any(s<=0)) {stop("s must be positive")}
  if (any(t0<0)) {stop("Non-decision time, t0, has to be non-negative")}
  if (any(st0<0)) {stop("Non-decision time range, st0, has to be non-negative")}

  if (!time_scaled) {
    wrt <- 0
    wint <- 0
    wx <- 1
  }

  pars <- prepare_RaceModel_parameter(response = 1,
                                      mu1, mu2,
                                      a, b,
                                      s, 0, 1,
                                      t0, st0,
                                      wx, wrt, wint,
                                      n)
  res <- matrix(NA, nrow=n, ncol=4)
  for (i in seq_len(length(pars$parameter_indices))) {
    ok_rows <- pars$parameter_indices[[i]]
    current_n <- length(ok_rows)
    out <- r_RM(current_n, pars$params[ok_rows[1],1:5], indep=TRUE, delta = delta, maxT = maxrt)

    ws <- pars$params[ok_rows[1],9:11]
    ws <- ws/sum(ws)

    # Since the Cpp function uses a different parametrization, -out[,3] is
    # exactly the distance of the loosing accumulator from its boundary
    if (time_scaled) {
      res[ok_rows,4] <- -ws[1]*out[,3] + ws[2]/sqrt(out[,1]) - ws[3]*out[,3]/sqrt(out[,1])
    } else {
      res[ok_rows,4] <- -out[,3]
    }
    # Since the Cpp function uses a different parametrization, we have to add
    # a (or b), i.e. -params[3] or -params[4] to the output to get the right
    # "state"
    out[,3] <- out[,3] - ifelse(out[,2]==1, pars$params[ok_rows[1],4], pars$params[ok_rows[1],3])
    ## Add the non-decision time to the rt-outcome
    out[,1] <- out[,1] + pars$params[ok_rows[1],12] + runif(current_n, 0, pars$params[ok_rows[1],8])
    res[ok_rows,1:3] <- out

  }

  res <- as.data.frame(res)
  names(res) <- c("rt", "response", "xl", "conf")

  return(res)
}

#' @rdname RaceModels
#' @export
rPCRM <- function (n, mu1, mu2, a, b,
                  wx=1, wrt=0, wint=0,
                  t0=0, st0=0, s=1,
                  time_scaled = TRUE, step_width=NULL,
                  delta=0.01, maxrt=15)
{
  if (any(missing(mu1), missing(mu2),
          missing(a), missing(b))) stop("mu1, mu2, a, and b must be supplied")
  if (any(c(a<=0, b<=0))) {stop("Both thresholds (a  and b) must be positive")}
  if (any(s<=0)) {stop("s must be positive")}
  if (any(t0<0)) {stop("Non-decision time, t0, has to be non-negative")}
  if (any(st0<0)) {stop("Non-decision time range, st0, has to be non-negative")}

  if (!time_scaled) {
    wrt <- 0
    wint <- 0
    wx <- 1
  }

  pars <- prepare_RaceModel_parameter(response = 1,
                                      mu1, mu2,
                                      a, b,
                                      s, 0, 1,
                                      t0, st0,
                                      wx, wrt, wint,
                                      n)
  res <- matrix(NA, nrow=n, ncol=4)
  for (i in seq_len(length(pars$parameter_indices))) {
    ok_rows <- pars$parameter_indices[[i]]
    current_n <- length(ok_rows)
    out <- r_RM(current_n, pars$params[ok_rows[1],1:5], indep=FALSE, delta = delta, maxT = maxrt)

    ws <- pars$params[ok_rows[1],9:11]
    ws <- ws/sum(ws)

    # Since the Cpp function uses a different parametrization, -out[,3] is
    # exactly the distance of the loosing accumulator from its boundary
    if (time_scaled) {
      res[ok_rows,4] <- -ws[1]*out[,3] + ws[2]/sqrt(out[,1]) - ws[3]*out[,3]/sqrt(out[,1])
    } else {
      res[ok_rows,4] <- -out[,3]
    }
    # Since the Cpp function uses a different parametrization, we have to add
    # a (or b), i.e. -params[3] or -params[4] to the output to get the right
    # "state"
    out[,3] <- out[,3] - ifelse(out[,2]==1, pars$params[ok_rows[1],4], pars$params[ok_rows[1],3])
    ## Add the non-decision time to the rt-outcome
    out[,1] <- out[,1] + pars$params[ok_rows[1],12] + runif(current_n, 0, pars$params[ok_rows[1],8])
    res[ok_rows,1:3] <- out

  }

  res <- as.data.frame(res)
  names(res) <- c("rt", "response", "xl", "conf")

  return(res)
}





prepare_RaceModel_parameter <- function(response,mu1, mu2,
                                        a, b, s, th1, th2,t0, st0, wx, wrt, wint, nn) {

  if ( (length(mu1) == 1) &
       (length(mu2) == 1) &
       (length(a) == 1) &
       (length(b) == 1) &
       (length(s) == 1) &
       (length(th1) == 1) &
       (length(th2) == 1) &
       (length(t0) == 1) &
       (length(st0) == 1) &
       (length(wx) == 1) &
       (length(wrt) == 1) &
       (length(wint) == 1)) {
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
    mu1 <- rep(mu1, length.out = nn)
    mu2 <- rep(mu2, length.out = nn)
    b <- rep(b, length.out = nn)
    s <- rep(s, length.out = nn)
    a <- rep(a, length.out = nn)
    th1 <- rep(th1, length.out = nn)
    th2 <- rep(th2, length.out = nn)
    t0 <- rep(t0, length.out = nn)
    st0 <- rep(st0, length.out = nn)
    wx <- rep(wx, length.out = nn)
    wrt <- rep(wrt, length.out = nn)
    wint <- rep(wint, length.out = nn)
  }
  th1[th1==-Inf] <- 0
  th2[th2==Inf] <- .Machine$double.xmax
  # Build parameter matrix (and divide a, v, and sv, by s)
  params <- cbind (mu1, mu2, -a, -b, s, th1, th2, st0, wx, wrt, wint, t0, numeric_bounds)

  # Check for illegal parameter values
  if(ncol(params)<13) stop("Not enough parameters supplied: probable attempt to pass NULL values?")
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

