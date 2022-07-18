#' Drift Diffusion Model with time-dependent confidence
#'
#' Likelihood function and random number generator for the Drift Diffuson Model
#' with confidence computet as 1/sqrt(DecisionTime). It includes following parameters:
#' DDM parameters: \code{a} (threshold separation), \code{z}
#' (starting point; relative), \code{v} (drift rate), \code{t0} (non-decision time/
#' response time constant), \code{d} (differences in speed of response execution),
#' \code{sv} (inter-trial-variability of drift), \code{st0} (inter-trial-variability
#' of non-decisional components), \code{sz} (inter-trial-variability of relative
#' starting point), \code{s} (diffusion constant).
#'
#' For the confidence part: \code{th1} and \code{th2} (lower and upper
#' thresholds for confidence interval).
#'
#' \strong{Note that the parameterization or defaults of non-decision time variability
#' \code{st0} and diffusion constant \code{s} differ from what is often found in the
#' literature.}
#'
#' @param rt a vector of RTs. Or for convenience also a \code{data.frame} with columns
#' \code{rt} and \code{response}.
#' @param response character vector, indicating the decision, i.e. which boundary was
#' met first. Possible values are \code{c("upper", "lower")} (possibly abbreviated) and
#' \code{"upper"} being the default. Alternatively, a numeric vector with values 1=lower
#' and 2=upper or -1=lower and 1=upper, respectively. For convenience, \code{response} is
#' converted via \code{as.numeric} also
#' allowing factors. Ignored if the first argument is a \code{data.frame}.
#' @param th1 together with th2: scalars or numerical vectors giving the lower and upper
#' bound of the interval, in which the accumulator should end at the time of the
#' confidence judgment (i.e. at time \code{rt}+\code{tau}). Only values with
#' \code{th2}>=\code{th1} are accepted.
#' @param th2 (see th1)
#'
#' @param a threshold separation. Amount of information that is considered for a decision.
#' Large values indicate a conservative decisional style. Typical range: 0.5 < \code{a} < 2
#' @param v drift rate. Average slope of the information accumulation process. The drift
#' gives information about the speed and direction of the accumulation of information.
#' Large (absolute) values of drift indicate a good performance. If received information
#' supports the response linked to the upper threshold the sign will be positive and vice
#' versa. Typical range: -5 < \code{v} < 5
#' @param t0 non-decision time or response time constant (in seconds). Lower bound for the
#' duration of all non-decisional processes (encoding and response execution). Typical
#' range: 0.1 < \code{t0} < 0.5. Default is 0.
#' @param z (by default relative) starting point. Indicator of an a priori bias in decision
#' making. When the relative starting point \code{z} deviates from \code{0.5}, the amount
#' of information necessary for a decision differs between response alternatives. Default
#' is \code{0.5} (i.e., no bias).
#' @param d differences in speed of response execution (in seconds). Positive values
#' indicate that response execution is faster for responses linked to the upper threshold
#' than for responses linked to the lower threshold. Typical range: -0.1 < \code{d} < 0.1.
#' Default is 0.
#' @param sz inter-trial-variability of starting point. Range of a uniform distribution
#' with mean \code{z} describing the distribution of actual starting points from specific
#' trials. Values different from 0 can predict fast errors (but can slow computation
#' considerably). Typical range: 0 < \code{sz} < 0.2. Default is 0. (Given in relative
#' range i.e. bounded by 2*min(z, 1-z))
#' @param sv inter-trial-variability of drift rate. Standard deviation of a normal
#' distribution with mean \code{v} describing the distribution of actual drift rates
#' from specific trials. Values different from 0 can predict slow errors.
#' Typical range: 0 < \code{sv} < 2. Default is 0.
#' @param st0 inter-trial-variability of non-decisional components. Range of a uniform
#' distribution with mean \code{t0 + st0/2} describing the distribution of actual
#' \code{t0} values across trials. Accounts for response times below \code{t0}.
#' Reduces skew of predicted RT distributions. Values different from 0 can slow computation
#' considerably. Typical range: 0 < \code{st0} < 0.2. Default is 0.
#' @param s diffusion constant. Standard deviation of the random noise of the diffusion
#' process (i.e., within-trial variability), scales \code{a}, \code{v}, \code{sv},
#' and \code{th}'s. Needs to be fixed to a constant in most applications. Default is 1.
#' Note that the default used by Ratcliff and in other applications is often 0.1.
#'
#' @param precision \code{numerical} scalar value. Precision of calculation. Corresponds
#' to the stepsize of integration w.r.t. z and t0. Default is 1e-5.
#' @param z_absolute logical. Determines whether z is treated as absolute start point
#' (TRUE) or relative (FALSE; default) to a.
#' @param stop_on_error Should the diffusion functions return 0 if the parameters values
#' are outside the allowed range (= \code{FALSE}) or produce an error in this case
#' (= \code{TRUE}).
#' @param stop_on_zero Should the computation of densities stop as soon as a density value of 0 occurs.
#' This may save a lot of time if the function is used for a likelihood function. Default: FALSE
#'
#' @param n integer. The number of samples generated.
#' @param delta numeric. Discretization step size for simulations in the stochastic process
#' @param maxrt numeric. Maximum decision time returned. If the simulation of the stochastic
#' process exceeds a decision time of `maxrt`, the `response` will be set to 0 and the `maxrt`
#' will be returned as `rt`.
#'
#' @return \code{dDDMConf} gives the density/likelihood/probability of the diffusion process
#' producing a decision of \code{response} at time \code{rt} and a confidence
#' judgment corresponding to the interval \[ \code{th1}, \code{th2}\].
#' The value will be a numeric vector of the same length as \code{rt}.
#'
#' \code{rDDMConf} returns a `data.frame` with three columns and `n` rows. Column names are `rt` (response
#' time), `response` (-1 (lower) or 1 (upper), indicating which bound was hit), and `conf` (the
#' value of the confidence measure; not discretized!).
#'
#' The distribution parameters (as well as \code{response}, \code{th1}
#' and \code{th2}) are recycled to the length of the result. In other words, the functions
#' are completely vectorized for all parameters and even the response boundary.
#'
#' @details The Ratcliff diffusion model (Ratcliff
#' and McKoon, 2008) is a mathematical model for two-choice discrimination tasks. It is
#' based on the assumption that information is accumulated continuously until one of two
#' decision thresholds is hit. For introduction see Ratcliff and McKoon (2008).
#'
#' This model incorporates the idea, that the decision time is informative for
#' stimulus difficulty and thus confidence is computed as a monotone function
#' of 1/sqrt(DecisionTime). Here, we use a given interval, given by \code{th1}
#' and \code{th2}, assuming that the data is given with discrete judgments and
#' preprocessed, s.t. these discrete ratings are translated to the respective intervals.
#'
#' All functions are fully vectorized across all parameters as well as the response to
#' match the length or \code{rt} (i.e., the output is always of length equal to \code{rt}).
#' This allows for trial wise parameters for each model parameter.
#'
#' For convenience, the function allows that the first argument is a \code{data.frame}
#' containing the information of the first and second argument in two columns (i.e.,
#' \code{rt} and \code{response}). Other columns (as well as passing \code{response}
#' separately argument) will be ignored.
#'
#' @note The parameterization of the non-decisional components, \code{t0} and \code{st0},
#' differs from the parameterization sometimes used in the literature.
#' In the present case \code{t0} is the lower bound of the uniform distribution of length
#' \code{st0}, but \emph{not} its midpoint. The parameterization employed here is in line
#' with the functions in the `rtdists` package.
#'
#' The default diffusion constant \code{s} is 1 and not 0.1 as in most applications of
#' Roger Ratcliff and others. Usually \code{s} is not specified as the other parameters:
#' \code{a}, \code{v}, and \code{sv}, may be scaled to produce the same distributions
#' (as is done in the code).
#'
#' The function code is basically an extension of the \code{ddiffusion} function from the
#' package \code{rtdists} for the Ratcliff diffusion model.
#'
#' @references Ratcliff, R., & McKoon, G. (2008). The diffusion decision model: Theory and data for two-choice decision tasks. \emph{Neural Computation}, 20(4), 873-922.
#'
#'
#' @author For the original rtdists package: Underlying C code by Jochen Voss and Andreas Voss. Porting and R wrapping by Matthew Gretton, Andrew Heathcote, Scott Brown, and Henrik Singmann. \code{qdiffusion} by Henrik Singmann. For the dDDMConf function the C code was extended by Sebastian Hellmann.
#'
#' @useDynLib dynConfiR, .registration = TRUE
#'
#' @name dDDMConf
#' @aliases DDMConf rDDMConf
#' @importFrom Rcpp evalCpp
#'
#' @examples
#' # Plot rt distribution ignoring confidence
#' curve(dDDMConf(x, "upper", 0, Inf, a=2, v=0.4, sz=0.2, sv=0.9), xlim=c(0, 2), lty=2)
#' curve(dDDMConf(x, "lower", 0, Inf, a=2, v=0.4, sz=0.2, sv=0.9), col="red", lty=2, add=TRUE)
#' curve(dDDMConf(x, "upper", 0, Inf, a=2, v=0.4),add=TRUE)
#' curve(dDDMConf(x, "lower", 0, Inf, a=2, v=0.4), col="red", add=TRUE)
#' # Generate a random sample
#' dfu <- rDDMConf(5000, a=2,v=0.5,t0=0,z=0.5,d=0,sz=0,sv=0, st0=2, s=1)
#' # Same RT distribution but upper and lower responses changed
#' dfl <- rDDMConf(50, a=2,v=-0.5,t0=0,z=0.5,d=0,sz=0,sv=0, st0=2, s=1)
#' head(dfu)
#'
#' dDDMConf(dfu, th1=0.5, th2=2.5, a=2, v=.5, st0=2)[1:5]
#' # Scaling diffusion parameters leads do same density values
#' s <- 2
#' dDDMConf(dfu, th1=0.5, th2=2.5, a=2*s, v=.5*s, s=2, st0=2)[1:5]
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   require(ggplot2)
#'   ggplot(dfu, aes(x=rt, y=conf))+
#'     stat_density_2d(aes(fill = ..density..), geom = "raster", contour = FALSE) +
#'     facet_wrap(~response)
#' }
#' boxplot(conf~response, data=dfu)
#'
#' # Restricting to specific confidence region
#' dfu <- dfu[dfu$conf >0 & dfu$conf <1,]
#' dDDMConf(dfu, th1=0, th2=1, a=2, v=0.5)[1:5]
#'
#' # If lower confidence threshold is higher than the upper, the function throws an error,
#' # except when stop_on_error is FALSE
#' dDDMConf(dfu[1:5,], th1=1, th2=0, a=2, v=0.5, stop_on_error = FALSE)
#'




#' @rdname dDDMConf
#' @export
dDDMConf <- function (rt, response="upper", th1,th2,tau=1,a,v,t0=0,z=0.5,d=0,sz=0,sv=0, st0=1,s=1,
                   precision=1e-5, z_absolute = FALSE,
                   stop_on_error=TRUE, stop_on_zero = FALSE)
{
  # for convenience accept data.frame as first argument.
  if (is.data.frame(rt)) {
    response <- rt$response
    rt <- rt$rt
  }

  nn <- length(rt)
  if (s != 1) {
    a <- a/s
    v <- v/s
    sv <- sv/s
    s <- 1
  }
  pars <- prepare_2DSD_parameter(response = response,
                                 a = a, v = v, t0 = t0, z = z,
                                 d = d, sz = sz, sv = sv, st0 = st0,
                                 tau=1, th1=th1, th2=th2,
                                 s = s, nn = nn, z_absolute = z_absolute,
                                 stop_on_error = stop_on_error)

  densities <- vector("numeric",length=nn)
  for (i in seq_len(length(pars$parameter_indices))) {
    ok_rows <- pars$parameter_indices[[i]]

    densities[ok_rows] <- d_DDMConf (rt[ok_rows],
                                  pars$params[ok_rows[1],1:11],
                                  precision,
                                  pars$params[ok_rows[1],12],
                                  stop_on_error, as.numeric(stop_on_zero))
  }
  abs(densities)
}

#' @rdname dDDMConf
#' @export
rDDMConf <- function (n, a,v,t0=0,z=0.5,d=0,sz=0,sv=0, st0=2,
                    s=1, delta=0.01, maxrt=15,
                   z_absolute = FALSE,  stop_on_error=TRUE)
{
  if (any(missing(a), missing(v))) stop("a and v must be supplied")

  pars <- prepare_2DSD_parameter(response = 1L,
                                 a = a, v = v, t0 = t0, z = z,
                                 d = d, sz = sz, sv = sv, st0 = st0,
                                 tau=1, th1=0, th2=1,
                                 s = s, nn = n, z_absolute = z_absolute,
                                 stop_on_error = stop_on_error)
  res <- matrix(NA, nrow=n, ncol=3)
  for (i in seq_len(length(pars$parameter_indices))) {
    ok_rows <- pars$parameter_indices[[i]]
    current_n <- length(ok_rows)
    out <- r_DDMConf(current_n, pars$params[ok_rows[1], 1:8],
                delta = delta, maxT =maxrt, stop_on_error)
    res[ok_rows,] <- out
  }

  res <- as.data.frame(res)
  names(res) <- c("rt", "response", "conf")
  return(res)
}



## Function to dynamical unload the .dll file
.onUnload <- function (libpath) {
  library.dynam.unload("dynConfiR", libpath)
}

