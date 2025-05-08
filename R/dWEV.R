#' Dynamical visibility, time, and evidence model (dynaViTE) and Dynamical weighted evidence and visibility model (dynWEV)
#'
#' Likelihood function and random number generator for the dynaViTE and dynWEV model
#' (Hellmann et al., 2023).
#' It includes following parameters from the drift diffusion model:
#' \code{a} (threshold separation),
#' \code{z} (starting point; relative),
#' \code{v} (drift rate),
#' \code{t0} (non-decision time/response time constant),
#' \code{d} (differences in speed of response execution),
#' \code{sv} (inter-trial-variability of drift),
#' \code{st0} (inter-trial-variability of non-decisional components),
#' \code{sz} (inter-trial-variability of relative starting point) and
#' \code{s} (diffusion constant).
#' For the computation of confidence following parameters were added:
#' \code{tau} (post-decisional accumulation time),
#' \code{w} (weight on the decision evidence (weight on visibility is (1-w))),
#' \code{muvis} (mean drift rate of visibility process),
#' \code{svis} (diffusion constant of visibility process),
#' \code{sigvis} (variability in drift rate of visibility accumulator),
#' \code{th1} and \code{th2} (lower and upper thresholds for confidence interval).
#' \code{lambda} for dynaViTE only, the exponent of judgment time for the division by judgment time in the confidence measure, and
#'  \strong{Note that the parametrization or defaults of non-decision time variability
#'  \code{st0} and diffusion constant \code{s} differ from what is often found in the literature.}
#'
#' @param rt a vector of RTs. Or for convenience also a \code{data.frame} with columns \code{rt}
#' and \code{response}.
#' @param response character vector, indicating the decision, i.e. which boundary was
#' met first. Possible values are \code{c("upper", "lower")} (possibly abbreviated) and
#' \code{"upper"} being the default. Alternatively, a numeric vector with values 1=lower
#' and 2=upper or -1=lower and 1=upper, respectively. For convenience, \code{response} is
#' converted via \code{as.numeric} also
#' allowing factors. Ignored if the first argument is a \code{data.frame}.
#' @param th1 together with `th2`: scalars or numerical vectors giving the lower and upper bound of the
#' interval of the confidence measure (see Details). Only values with \code{th2}>=\code{th1} are accepted.
#' @param th2 (see `th1`)
#'
#' @param a threshold separation. Amount of information that is considered for a decision. Large values
#' indicate a conservative decisional style. Typical range: 0.5 < \code{a} < 2
#' @param v drift rate of decision process. Average slope of the information accumulation
#' process. The drift gives information about the speed and direction of the accumulation of information.
#' Large (absolute) values of drift indicate a good performance. If received information supports the
#' response linked to the upper threshold the sign will be positive and vice versa.
#' Typical range: -5 < \code{v} < 5
#' @param t0 non-decision time or response time constant (in seconds). Lower bound for the duration of all
#' non-decisional processes (encoding and response execution).
#' Typical range: 0.1 < \code{t0} < 0.5. Default is 0.
#' @param z (by default relative) starting point of decision process. Indicator of an a priori bias in decision
#' making. When the relative starting point \code{z} deviates from \code{0.5}, the amount
#' of information necessary for a decision differs between response alternatives. Default
#' is \code{0.5} (i.e., no bias).
#' @param d differences in speed of response execution (in seconds). Positive values indicate that
#' response execution is faster for responses linked to the upper threshold than for responses linked to
#' the lower threshold. Typical range: -0.1 < \code{d} < 0.1. Default is 0.
#' @param sz inter-trial-variability of starting point. Range of a uniform distribution with mean
#' \code{z} describing the distribution of actual starting points from specific trials. Values different
#' from 0 can predict fast errors (but can slow computation considerably). Typical range:
#' 0 < \code{sz} < 0.2. Default is 0. (Given in relative range i.e. bounded by 2*min(z, 1-z))
#' @param sv inter-trial-variability of drift rate of decision process. Standard deviation of a normal
#' distribution with mean \code{v} describing the distribution of actual drift rates from specific trials.
#' Values different from 0 can predict slow errors. Typical range: 0 < \code{sv} < 2. Default is 0.
#' @param st0 inter-trial-variability of non-decisional components. Range of a uniform distribution with
#' mean \code{t0 + st0/2} describing the distribution of actual \code{t0} values across trials.
#' Accounts for response times below \code{t0}. Reduces skew of predicted RT distributions.
#' Values different from 0 can slow computation considerably. Typical range: 0 < \code{st0} < 0.2.
#' Default is 0.
#' @param s diffusion constant of decision process; standard deviation of the random noise of the
#' diffusion process (i.e., within-trial variability), scales other parameters (see Note). Needs to
#' be fixed to a constant in most applications. Default is 1. Note that the default used by Ratcliff
#' and in other applications is often 0.1.
#' @param tau post-decisional accumulation time; the length of the time period after the decision was
#' made until the confidence judgment is made. Range: \code{tau}>0. Default: \code{tau}=1.
#' @param w weight put on the final state of the decision accumulator for confidence computation.
#' 1-w is the weight on the visibility accumulator. Range: 0<\code{w}<1. Default: \code{w}=0.5.
#' @param muvis mean drift of visibility process. If `NULL` (default), `muvis` will be set to the
#' absolute value of `v`.
#' @param svis diffusion constant of visibility process. Range: \code{svis}>0. Default: \code{svis}=1.
#' @param sigvis the variability in drift rate of the visibility process (which varies independently
#' from the drift rate in decision process). Range: \code{sigvis}>=0. Default: \code{sigvis}=0.
#' @param lambda power for judgment time in the division of the confidence measure
#' by the judgment time (Default: 0, i.e. no division which is the version of
#' dynWEV proposed by Hellmann et al., 2023)
#'
#' @param  simult_conf logical. Whether in the experiment confidence was reported simultaneously
#' with the decision. If that is the case decision and confidence judgment are assumed to have happened
#' subsequent before the response. Therefore `tau` is included in the response time. If the decision was
#' reported before the confidence report, `simul_conf` should be `FALSE`.
#' @param precision numerical scalar value. Precision of calculation. Determines the
#' step size of integration w.r.t. `z` and `t0`.
#' Represents the number of decimals precisely computed on average. Default is 6.
#' @param z_absolute logical. Determines whether `z` is treated as absolute start point
#' (`TRUE`) or relative (`FALSE`; default) to `a`.
#' @param stop_on_error Should the diffusion functions return 0 if the parameters values are
#' outside the allowed range (= \code{FALSE}) or produce an error in this case (= \code{TRUE}).
#' @param stop_on_zero Should the computation of densities stop as soon as a density value of 0 occurs.
#' This may save a lot of time if the function is used for a likelihood function. Default: FALSE
#' @param n integer. The number of samples generated.
#' @param delta numeric. Discretization step size for simulations in the stochastic process
#' @param maxrt numeric. Maximum decision time returned. If the simulation of the stochastic
#' process exceeds a decision time of `maxrt`, the `response` will be set to 0 and the `maxrt`
#' will be returned as `rt`.
#' @param process_results logical. Whether the output simulations should contain the final
#' state of the decision (and visibility) process as additional column. Default is FALSE, meaning that
#' no additional columns for the final process states are returned.
#'
#' @return \code{ddynaViTE} gives the density/likelihood/probability of the diffusion process producing
#' a decision of \code{response} at time \code{rt} and a confidence judgment corresponding to the
#' interval \[ \code{th1}, \code{th2}\]. The value will be a numeric vector of the same length as
#' \code{rt}.
#'
#' \code{rdynaViTE} returns a `data.frame` with three columns and `n` rows. Column names are `rt` (response
#' time), `response` (-1 (lower) or 1 (upper), indicating which bound was hit), and `conf` (the
#' value of the confidence measure; not discretized!).
#'
#' The distribution parameters (as well as \code{response}, \code{tau}, \code{th1} and \code{th2},
#' \code{w} and \code{sig}) are recycled to the length of the result. In other words, the functions
#' are completely vectorized for all parameters and even the response boundary.
#'
#' @details
#' The function dWEV was renamed to ddynaViTE in version 0.1.0 of the package.
#' It is still here for reasons of backwards compatibility. The function just calls the
#' ddynaViTE function (and produces a deprecation warning).
#'
#' The dynamical visibility, time, and evidence (dynaViTE) model
#' and the weighted evidence and visibility model (dynWEV) are extensions of the 2DSD model
#' for decision confidence (see \code{\link{d2DSD}}). It assumes that the decision follows a drift
#' diffusion model with two additional assumptions to account for confidence. First, there is a
#' post-decisional period of further evidence accumulation `tau`. Second, another accumulation process
#' accrues information about stimulus reliability (the visibility process) including also evidence
#' about decision irrelevant features. See Hellmann et al. (2023) for more information.
#' The measure for confidence is then a weighted sum of the final state of the decision process X
#' and the visibility process V over a power-function of total accumulation time,
#' i.e. for a decision time T (which is not the response time), the confidence variable is
#' \deqn{conf = \frac{wX(T+\tau) + (1-w) V(T+\tau)}{(T+\tau)^\lambda}.}
#' The dynWEV model is a special case of dynaViTE, with the parameter `lambda=0`.
#'
#' All functions are fully vectorized across all parameters as well as the response to match the
#' length or \code{rt} (i.e., the output is always of length equal to \code{rt}). This allows for
#' trial wise parameters for each model parameter.
#'
#' For convenience, the function allows that the first argument is a \code{data.frame} containing
#' the information of the first and second argument in two columns (i.e., \code{rt} and \code{response}).
#' Other columns (as well as passing \code{response} separately argument) will be ignored.
#'
#' @note The parameterization of the non-decisional components, \code{t0} and \code{st0},
#' differs from the parameterization sometimes used in the literature.
#' In the present case \code{t0} is the lower bound of the uniform distribution of length
#' \code{st0}, but \emph{not} its midpoint. The parameterization employed here is in line
#' with the functions in the `rtdists` package.
#'
#' The default diffusion constant \code{s} is 1 and not 0.1 as in most applications of
#' Roger Ratcliff and others. Usually \code{s} is not specified as the other parameters:
#' \code{a}, \code{v}, \code{sv}, \code{muvis}, \code{sigvis}, and \code{svis} respectively,
#' may be scaled to produce the same distributions (as is done in the code).
#'
#' The function code is basically an extension of the \code{ddiffusion} function from the
#' package \code{rtdists} for the Ratcliff diffusion model.
#'
#' @references  Hellmann, S., Zehetleitner, M., & Rausch, M. (2023). Simultaneous modeling of choice, confidence and response time in visual perception. \emph{Psychological Review} 2023 Mar 13. doi: 10.1037/rev0000411. Epub ahead of print. PMID: 36913292.
#'
#' @author Sebastian Hellmann
#'
#' @useDynLib dynConfiR, .registration = TRUE
#'
#' @name dynaViTE
#' @aliases WEVmodel ddynaViTE ddynWEV rdynaViTE dynWEV dWEV
#' @importFrom Rcpp evalCpp
#'
#' @examples
#' # Plot rt distribution ignoring confidence
#' curve(ddynaViTE(x, "upper", -Inf, Inf, tau=1, a=2, v=0.4, sz=0.2, sv=0.9), xlim=c(0, 2), lty=2)
#' curve(ddynaViTE(x, "lower", -Inf, Inf,tau=1, a=2, v=0.4, sz=0.2, sv=0.9))
#' curve(ddynaViTE(x, "upper", -Inf, Inf,  tau=1, a=2, v=0.4),add=TRUE)
#' curve(ddynaViTE(x, "lower", -Inf, Inf, tau=1, a=2, v=0.4), col="red", add=TRUE)
#' # Generate a random sample
#' df1 <- rdynaViTE(5000, a=2,v=0.5,t0=0,z=0.5,d=0,sz=0,sv=0, st0=0,  tau=1, s=1, w=0.9)
#' # Same RT and response distribution but different confidence distribution
#' df2 <- rdynaViTE(5000, a=2,v=0.5,t0=0,z=0.5,d=0,sz=0,sv=0, st0=0,  tau=1, s=1, w=0.1)
#' head(df1)
#'
#' # Scaling diffusion parameters leads do same density values
#' ddynaViTE(df1[1:5,], th1=-Inf, th2=Inf, a=2, v=.5)[1:5]
#' s <- 2
#' ddynaViTE(df1[1:5,], th1=-Inf, th2=Inf, a=2*s, v=.5*s, s=2)[1:5]
#'
#' # Diffusion constant also scales confidence parameters
#' ddynaViTE(df1[1:5,], th1=0.2, th2=1, a=2, v=.5, sv=0.2, w=0.5, sigvis = 0.2, svis = 1)[1:5]
#' s <- 2
#' ddynaViTE(df1[1:5,], th1=0.2*s, th2=1*s, a=2*s, v=.5*s, s=2,
#'      sv=0.2*s, w=0.5, sigvis=0.2*s, svis=1*s)[1:5]
#'
#'
#' two_samples <- rbind(cbind(df1, w="high"),
#'                      cbind(df2, w="low"))
#' # no difference in RT distributions
#' boxplot(rt~w+response, data=two_samples)
#' # but different confidence distributions
#' boxplot(conf~w+response, data=two_samples)
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   require(ggplot2)
#'   ggplot(two_samples, aes(x=rt, y=conf))+
#'     stat_density_2d(aes(fill = after_stat(density)), geom = "raster", contour = FALSE) +
#'     xlim(c(0, 2))+ ylim(c(-1.5, 4))+
#'     facet_grid(cols=vars(w), rows=vars(response), labeller = "label_both")
#' }
#'
#' # Restricting to specific confidence region
#' df1 <- df1[df1$conf >0 & df1$conf <1,]
#' ddynaViTE(df1[1:5,], th1=0, th2=1, a=2, v=0.5)[1:5]
#'
#' # If lower confidence threshold is higher than the upper, the function throws an error,
#' # except when stop_on_error is FALSE
#' ddynaViTE(df1[1:5,], th1=1, th2=0, a=2, v=0.5, stop_on_error = FALSE)
#'

#' @rdname dynaViTE
#' @export
dWEV <- function (rt, response="upper", th1,th2, a,v,t0=0,z=0.5,d=0,sz=0,sv=0, st0=0,
                       tau=1, w=0.5, muvis=NULL, sigvis=0, svis=1,
                       lambda = 0, s=1,  simult_conf = FALSE, precision=6, z_absolute = FALSE,
                       stop_on_error=TRUE, stop_on_zero=FALSE) {
  warning("The function dWEV was renamed to ddynaViTE in version 0.1.0 of the package.\\
          It will be removed in future releases of the package!")
  return(ddynaViTE(rt, response="upper", th1,th2, a,v,t0=0,z=0.5,d=0,sz=0,sv=0, st0=0,
                         tau=1, w=0.5, muvis=NULL, sigvis=0, svis=1,
                         lambda = 0, s=1,  simult_conf = FALSE, precision=6, z_absolute = FALSE,
                         stop_on_error=TRUE, stop_on_zero=FALSE))
}

#' @rdname dynaViTE
#' @export
ddynaViTE <- function (rt, response="upper", th1,th2, a,v,t0=0,z=0.5,d=0,sz=0,sv=0, st0=0,
                  tau=1, w=0.5, muvis=NULL, sigvis=0, svis=1,
                  lambda = 0, s=1,  simult_conf = FALSE, precision=6, z_absolute = FALSE,
                  stop_on_error=TRUE, stop_on_zero=FALSE)
{
  # for convenience accept data.frame as first argument.
  if (is.data.frame(rt)) {
    response <- rt$response
    rt <- rt$rt
  }

  nn <- length(rt)
  if (is.null(muvis)) {
    muvis <- abs(v)
  }
  pars <- prepare_WEV_parameter(response = response,
                                a = a, v = v, t0 = t0, z = z,
                                d = d, sz = sz, sv = sv, st0 = st0,
                                tau=tau, th1=th1, th2=th2,
                                w=w, muvis=muvis, svis=svis,sigvis=sigvis,
                                lambda=lambda,
                                s = s, nn = nn, z_absolute = z_absolute,
                                stop_on_error = stop_on_error)

  densities <- vector("numeric",length=nn)
  for (i in seq_len(length(pars$parameter_indices))) {
    ok_rows <- pars$parameter_indices[[i]]
    if (simult_conf) {
      # if confidence and decision are given simultaneously, subtract tau from rt
      # because both processes are assumed to happen subsequentially and therefore
      # observable response time is the sum of decision time, interjudgment time
      # and non-decision component
      rt[ok_rows] <- rt[ok_rows]-pars$params[ok_rows[1],9]
    }
    densities[ok_rows] <- d_WEVmu (rt[ok_rows],
                                   pars$params[ok_rows[1],1:16],
                                   precision,
                                   pars$params[ok_rows[1],17],
                                   stop_on_error, as.numeric(stop_on_zero))
  }
  abs(densities)
}


#' @rdname dynaViTE
#' @export
rdynaViTE <- function (n, a,v,t0=0,z=0.5,d=0,sz=0,sv=0, st0=0,
                  tau=1, w=0.5, muvis=NULL, sigvis=0, svis=1,
                  lambda=0, s=1, delta=0.01, maxrt=15, simult_conf = FALSE,
                  z_absolute = FALSE,  stop_on_error=TRUE, process_results=FALSE)
{
  if (is.null(muvis)) muvis <- abs(v)
  if (any(missing(a), missing(v))) stop("a and v must be supplied")

  pars <- prepare_WEV_parameter(response = 1L,
                                a = a, v = v, t0 = t0, z = z,
                                d = d, sz = sz, sv = sv, st0 = st0,
                                tau=tau, th1=0, th2=1,
                                w=w, muvis=muvis, svis=svis,sigvis=sigvis,
                                lambda=lambda,
                                s = s, nn = n, z_absolute = z_absolute,
                                stop_on_error = stop_on_error)
  res <- matrix(NA, nrow=n, ncol=6)
  for (i in seq_len(length(pars$parameter_indices))) {
    ok_rows <- pars$parameter_indices[[i]]
    current_n <- length(ok_rows)
    out <- r_WEV(current_n, pars$params[ok_rows[1], c(1:9, 12:16)],
                  delta = delta, maxT =maxrt, stop_on_error)

    out[,3] <- out[,3]/pars$params[ok_rows[1], 11]  # multiply by s (diffusion constant)
    out[,4] <- out[,4]/pars$params[ok_rows[1], 11]  # multiply by s (diffusion constant)
    out[,5] <- out[,5]/pars$params[ok_rows[1], 11]  # multiply by s (diffusion constant)
    if (simult_conf) {
      # add tau to response times, if choice and confidence given simultaneously
      out[,1] <- out[,1] + pars$params[ok_rows[1], 9]
    }
    res[ok_rows,] <- out
  }

  res <- as.data.frame(res)

  names(res) <- c("rt", "response", "conf", "dec", "vis", "mu")
  if (!process_results) res <- res[,1:3]
  return(res)
}


recalc_t0 <- function (t0, st0) { t0 <- t0 + st0/2 }

prepare_WEV_parameter <- function(response,
                                  a, v, t0, z, d,
                                  sz, sv, st0, tau, th1, th2,
                                  w, muvis, svis, sigvis,
                                  lambda,
                                  s, nn,
                                  z_absolute = FALSE,
                                  stop_on_error) {
  if(any(missing(a), missing(v))) stop("a and v must be supplied")
  if ( (length(s) == 1) &
       (length(a) == 1) &
       (length(v) == 1) &
       (length(t0) == 1) &
       (length(z) == 1) &
       (length(d) == 1) &
       (length(sz) == 1) &
       (length(sv) == 1) &
       (length(st0) == 1)&
       (length(tau) ==1) &
       (length(th1) ==1) &
       (length(th2) ==1) &
       (length(muvis) ==1) &
       (length(svis) ==1) &
       (length(w) == 1) &
       (length(sigvis) == 1)&
       (length(lambda) == 1)) {
    skip_checks <- TRUE
  } else {
    skip_checks <- FALSE
  }

  # Build parameter matrix
  # Convert boundaries to numeric if necessary
  if (all(response %in% c(-1, 1))) {
    response <- ifelse(response==1, 2, 1)
  }
  if (is.character(response)) {
    response <- match.arg(response, choices=c("upper", "lower"),several.ok = TRUE)
    numeric_bounds <- ifelse(response == "upper", 2L, 1L)
  }
  else {
    response <- as.numeric(response)
    if (any(!(response %in% 1:2)))
      stop("response needs to be either 'upper', 'lower', or as.numeric(response) %in% 1:2!")
    numeric_bounds <- as.integer(response)
  }

  numeric_bounds <- rep(numeric_bounds, length.out = nn)
  if (!skip_checks) {
    # all parameters brought to length of rt
    s <- rep(s, length.out = nn)
    a <- rep(a, length.out = nn)
    v <- rep(v, length.out = nn)
    t0 <- rep(t0, length.out = nn)
    z <- rep(z, length.out = nn)
    d <- rep(d, length.out = nn)
    sz <- rep(sz, length.out = nn)
    sv <- rep(sv, length.out = nn)
    st0 <- rep(st0, length.out = nn)
    tau <- rep(tau, length.out = nn)
    th1 <- rep(th1, length.out = nn)
    th2 <- rep(th2, length.out = nn)
    w <- rep(w, length.out = nn)
    svis <- rep(svis, length.out = nn)
    muvis <- rep(muvis, length.out = nn)
    sigvis <- rep(sigvis, length.out = nn)
    lambda <- rep(lambda, length.out = nn)
  }
  th1[th1==-Inf] <- - .Machine$double.xmax
  th2[th2==Inf] <- .Machine$double.xmax

  if (z_absolute) {
    z <- z/a  # transform z from absolute to relative scale (which is currently required by the C code)
    sz <- sz/a # transform sz from absolute to relative scale (which is currently required by the C code)
  }
  t0 <- recalc_t0 (t0, st0)

  # Build parameter matrix (and divide a, v, sv, sigvis and by s )

  params <- cbind (a/s, v/s, t0, d, sz, sv/s, st0, z,
                   tau, th1/s, th2/s, lambda, w, muvis/s, sigvis/s, svis/s, numeric_bounds)

  # Check for illegal parameter values
  if(ncol(params)<17) stop("Not enough parameters supplied: probable attempt to pass NULL values?")
  if(!is.numeric(params)) stop("Parameters need to be numeric.")
  if (any(is.na(params)) || !all(is.finite(params))) {
    if (stop_on_error) stop("Parameters need to be numeric and finite.")
    else return(rep(0, nn))
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


## Function to dynamical unload the .dll file
.onUnload <- function (libpath) {
  library.dynam.unload("dynConfiR", libpath)
}


