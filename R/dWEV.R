#' Dynamical weighted evidence and visibility model (dynWEV)
#'
#' Likelihood function for the dynWEV model (Hellmann & Rausch).
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
#' \code{tau} (postdecisional accumulation time),
#' \code{w} (weight on the decision evidence (weight on visibility is (1-w))),
#' \code{muvis} (mean drift rate of visibility process),
#' \code{svis} (diffusion constant of visibility process),
#' \code{sigvis} (variability in drift rate of visibility accumulator), and
#' \code{th1} and \code{th2} (lower and upper thresholds for confidence interval).
#'  \strong{Note that the parameterization or defaults of non-decision time variability
#'  \code{st0} and diffusion constant \code{s} differ from what is often found in the literature.}
#'  The code is basically an extension of the \code{\link[rtdists:Diffusion]{ddiffusion}}
#'  function for the Ratcliff diffusion model.
#'
#' @param rt a vector of RTs. Or for convenience also a \code{data.frame} with columns \code{rt} and \code{response} (such as returned from \code{rdiffusion} or \strong{rLBA}). See examples.
#' @param response character vector. Which response boundary should be tested? Possible values are \code{c("upper", "lower")}, possibly abbreviated and \code{"upper"} being the default. Alternatively, a numeric vector with values 1=lower and 2=upper. For convenience, \code{response} is converted via \code{as.numeric} also allowing factors (see examples). Ignored if the first argument is a \code{data.frame}.
#' @param th1 together with th2: scalars or numerical vectors giving the lower and upper bound of the interval, in which the accumulator should end at the time of the confidence judgement (i.e. at time \code{rt}+\code{tau}). Only values with \code{th2}>=\code{th1} are accepted.
#' @param th2 (see th1)
#'
#' @param a threshold separation. Amount of information that is considered for a decision. Large values indicate a conservative decisional style. Typical range: 0.5 < \code{a} < 2
#' @param v drift rate. Average slope of the information accumulation process. The drift gives information about the speed and direction of the accumulation of information. Large (absolute) values of drift indicate a good performance. If received information supports the response linked to the upper threshold the sign will be positive and vice versa. Typical range: -5 < \code{v} < 5
#' @param t0 non-decision time or response time constant (in seconds). Lower bound for the duration of all non-decisional processes (encoding and response execution). Typical range: 0.1 < \code{t0} < 0.5
#' @param z starting point. Indicator of an a priori bias in decision making. When the relative starting point \code{z} deviates from \code{0.5*a}, the amount of information necessary for a decision differs between response alternatives. Default is \code{0.5*a} (i.e., no bias).
#' @param d differences in speed of response execution (in seconds). Positive values indicate that response execution is faster for responses linked to the upper threshold than for responses linked to the lower threshold. Typical range: -0.1 < \code{d} < 0.1. Default is 0.
#' @param sz inter-trial-variability of starting point. Range of a uniform distribution with mean \code{z} describing the distribution of actual starting points from specific trials. Values different from 0 can predict fast errors (but can slow computation considerably). Typical range: 0 < \code{sz} < 0.2. Default is 0. (Given in relative range i.e. bounded by 2*min(z, 1-z))
#' @param sv inter-trial-variability of drift rate. Standard deviation of a normal distribution with mean \code{v} describing the distribution of actual drift rates from specific trials. Values different from 0 can predict slow errors. Typical range: 0 < \code{sv} < 2. Default is 0.
#' @param st0 inter-trial-variability of non-decisional components. Range of a uniform distribution with mean \code{t0 + st0/2} describing the distribution of actual \code{t0} values across trials. Accounts for response times below \code{t0}. Reduces skew of predicted RT distributions. Values different from 0 can slow computation considerably. Typical range: 0 < \code{st0} < 0.2. Default is 0.
#' @param s diffusion constant of decision process; standard deviation of the random noise of the diffusion process (i.e., within-trial variability), scales other parameters (see Note). Needs to be fixed to a constant in most applications. Default is 1. Note that the default used by Ratcliff and in other applications is often 0.1.
#' @param tau postdecisional accumulation time; the length of the time period after the decision was made until the confidence judgement is made. Range: \code{tau}>0. Default: \code{tau}=1.
#' @param w weight put on the final state of the decision accumulator. 1-w is the weight on the visibility accumulator. Range: 0<\code{w}<1. Default: \code{w}=0.5.
#' @param muvis mean drift of visibility accumulation process. Default: \code{muvis}=1.
#' @param svis diffusion constant of visibility accumulation process. Range: \code{svis}>0. Default: \code{svis}=1.
#' @param sigvis the variability in drift rate of the visibility accumulator (which varies independently from the drift rate in decision process). Range: \code{sigvis}>=0. Default: \code{sigvis}=0.
#'
#' @param  simult_conf logical. Whether in the experiment confidence was reported simultaneously
#' with the decision, as then decision and confidence judgment are assumed to have happened
#' subsequent before response and computations are different, when there is an observable
#' interjudgment time (then simult_conf should be FALSE).
#' @param precision \code{numerical} scalar value. Precision of calculation. #' @param precision \code{numerical} scalar value. Precision of calculation. Corresponds to the stepsize of integration w.r.t. z and t0. Default is 1e-5.
#' @param z_absolute logical. Determines whether z is treated as absolute start point (TRUE) or relative (FALSE; default) to a.
#' @param stop_on_error Should the diffusion functions return 0 if the parameters values are outside the allowed range (= \code{FALSE}) or produce an error in this case (= \code{TRUE}).
#'
#' @return \code{d2DSD(WEV(d/mu/absmu))} gives the density/likelihood/probability of the diffusion process producing a decision of \code{response} at time \code{rt} and a confidence judgment corresponding to the interval \[ \code{th1}, \code{th2}\]. The value will be a numeric vector of the same length as \code{rt}.
#'
#' The distribution parameters (as well as \code{response}, \code{tau}, \code{th1} and \code{th2}, \code{w} and \code{sig}) are recycled to the length of the result. In other words, the functions are completely vectorized for all parameters and even the response boundary.
#'
#' @details The function is a extension of the \code{ddiffusion} function from the package \code{rtdists} for the Ratcliff diffusion model. The Ratcliff diffusion model (Ratcliff, 1978) is a mathematical model for two-choice discrimination tasks. It is based on the assumption that information is accumulated continuously until one of two decision thresholds is hit. For introductions see Ratcliff and McKoon (2008), Voss, Rothermund, and Voss (2004), Voss, Nagler, and Lerche (2013), or Wagenmakers (2009).
#'
#' The dynWEV model is an extension of the drift diffusion model that incorporates the idea of an additional process accumulating information about the general visibility (also decision irrelevant features) of the stimulus as well as a postdecisional accumulation period (\strong{Rausch, Hellmann & Zehetleitner (2018)}).
#'
#' All functions are fully vectorized across all parameters as well as the response to match the length or \code{rt} (i.e., the output is always of length equal to \code{rt}). This allows for trialwise parameters for each model parameter.
#'
#' For convenience, the function allows that the first argument is a \code{data.frame} containing the information of the first and second argument in two columns (i.e., \code{rt} and \code{response}). Other columns (as well as passing \code{response} separately argument) will be ignored.
#'
#' @note The parameterization of the non-decisional components, \code{t0} and \code{st0}, differs from the parameterization used by, for example, Andreas Voss or Roger Ratcliff. In the present case \code{t0} is the lower bound of the uniform distribution of length \code{st0}, but \emph{not} its midpoint. The parameterization employed here is in line with the parametrization for the \strong{LBA} code (where \code{t0} is also the lower bound).
#'
#' The default diffusion constant \code{s} is 1 and not 0.1 as in most applications of Roger Ratcliff and others. Usually \code{s} is not specified as the other parameters:
#' \code{a}, \code{v}, \code{sv}, \code{muvis}, \code{sigvis}, and \code{svis} respectively, may be scaled to produce the same distributions (as is done in the code).
#'
#'
#' @references Ratcliff, R. (1978). A theory of memory retrieval. \emph{Psychological Review}, 85(2), 59-108.
#'
#' Ratcliff, R. (1978). A theory of memory retrieval. \emph{Psychological Review}, 85(2), 59-108.
#'
#' Ratcliff, R., & McKoon, G. (2008). The diffusion decision model: Theory and data for two-choice decision tasks. \emph{Neural Computation}, 20(4), 873-922.
#'
#' Rausch, M., Hellmann, S., & Zehetleitner, M. (2018). Confidence in masked orientation judgments is informed by both evidence and visibility. \emph{Attention, Perception, & Psychophysics}, 80(1), 134–154.  doi: 10.3758/s13414-017-1431-5
#'
#' Voss, A., Rothermund, K., & Voss, J. (2004). Interpreting the parameters of the diffusion model: An empirical validation. \emph{Memory & Cognition}. Vol 32(7), 32, 1206-1220.
#'
#' Voss, A., Nagler, M., & Lerche, V. (2013). Diffusion Models in Experimental Psychology: A Practical Introduction. \emph{Experimental Psychology}, 60(6), 385-402. doi:10.1027/1618-3169/a000218
#'
#' Wagenmakers, E.-J., van der Maas, H. L. J., & Grasman, R. P. P. P. (2007). An EZ-diffusion model for response time and accuracy. \emph{Psychonomic Bulletin & Review}, 14(1), 3-22.
#'
#' Wagenmakers, E.-J. (2009). Methodological and empirical developments for the Ratcliff diffusion model of response times and accuracy. \emph{European Journal of Cognitive Psychology}, 21(5), 641-671.
#'
#'
#' @author For the original rtdists package: Underlying C code by Jochen Voss and Andreas Voss. Porting and R wrapping by Matthew Gretton, Andrew Heathcote, Scott Brown, and Henrik Singmann. \code{qdiffusion} by Henrik Singmann. For the d2DSD function the C code was extended by Sebastian Hellmann.
#'
#' @useDynLib dynConfiR, .registration = TRUE
#'
#' @name dWEV
#' @aliases WEVmodel dynWEV
#' @importFrom Rcpp evalCpp
#'


# [MG 20150616]
# In line with LBA and rtdists, adjust t0 to be the lower bound of the non-decision time distribution rather than the average
# Called from prepare_WEV_parameter
recalc_t0 <- function (t0, st0) { t0 <- t0 + st0/2 }

prepare_WEV_parameter <- function(response,
                                  a, v, t0, z, d,
                                  sz, sv, st0, tau, th1, th2, w, muvis, svis, sigvis,
                                  s, nn,
                                  z_absolute = FALSE,
                                  stop_on_error) {
  if(any(missing(a), missing(v), missing(t0))) stop("a, v, and/or t0 must be supplied")
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
       (length(sigvis) == 1)) {
    skip_checks <- TRUE
  } else {
    skip_checks <- FALSE
  }

  # Build parameter matrix
  # Convert boundaries to numeric if necessary
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
                   tau, th1/s, th2/s, w, muvis/s, sigvis/s, svis/s, numeric_bounds)

  # Check for illegal parameter values
  if(ncol(params)<16) stop("Not enough parameters supplied: probable attempt to pass NULL values?")
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


#' @rdname dWEV
#' @export
dWEV <- function (rt,th1,th2,response="upper",a,v,t0,z=0.5,d=0,sz=0,sv=0, st0=0,
                  tau=1, w=0.5, muvis=NULL, sigvis=0, svis=1,
                  s=1, simult_conf = FALSE, precision=1e-5, z_absolute = FALSE,  stop_on_error=TRUE)
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
                                s = s, nn = nn, z_absolute = z_absolute,
                                stop_on_error = stop_on_error)

  densities <- vector("numeric",length=nn)
  for (i in seq_len(length(pars$parameter_indices))) {
    ok_rows <- pars$parameter_indices[[i]]
    if (simult_conf) {
      # if confidence and decision are given simultaneously, subtract tau from rt
      # because both processes are assumed to happen subsequentially and therefore
      # observable response time is the sum of decision time, inter-judgment time
      # and non-decision component
      rt[ok_rows] <- rt[ok_rows]-pars$params[ok_rows[1],9]
    }
    densities[ok_rows] <- d_WEVmu (rt[ok_rows],
                                   pars$params[ok_rows[1],1:15],
                                   precision,
                                   pars$params[ok_rows[1],16],
                                   stop_on_error)
  }
  abs(densities)
}

## Function to dynamical unload the .dll file
.onUnload <- function (libpath) {
  library.dynam.unload("dynConfiR", libpath)
}


