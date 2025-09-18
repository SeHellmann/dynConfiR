#' Simulation of confidence ratings and RTs in leaky competing accumulator model
#'
#' Simulates the decision responses, reaction times and state of the loosing accumulator
#' together with a confidence measure in the leaky competing accumulator model.
#' Optionally, there is a post-decisional accumulation period, where the processes continues.

#' @param n integer. number of samples.
#' @param mu1 mean momentary evidence for alternative 1
#' @param mu2 mean momentary evidence for alternative 2
#' @param th1 decision threshold for alternative 1
#' @param th2 decision threshold for alternative 2
#' @param k leakage (default: 0)
#' @param beta inhibition  (default: 0)
#' @param SPV variation in starting points  (default: 0)
#' @param tau fixed post decisional accumulation period  (default: 0)
#' @param wx weight on balance of evidence in confidence measure  (default: 1)
#' @param wrt weight on RT in confidence measure  (default: 0)
#' @param wint weight on interaction of evidence and RT in confidence measure (default: 0)
#' @param t0 minimal non-decision time (default: 0)
#' @param st0 range of uniform distribution of non-decision time (default: 0)
#' @param pi factor for input dependent noise of infinitesimal variance of processes (default: 0)
#' @param sig input independent component of infinitesimal variance of processes (default: 1)
#'
#' @param time_scaled logical. Whether a time_scaled transformation for the confidence measure should
#' be used.
#' @param  simult_conf logical. Whether in the experiment confidence was reported simultaneously
#' with the decision. If that is the case decision and confidence judgment are assumed to have happened
#' subsequent before the response. Therefore `tau` is included in the response time. If the decision was
#' reported before the confidence report, `simul_conf` should be `FALSE`.
#' @param delta numerical. Size of steps for the discretized simulation (see details).
#' @param maxrt numerical. Maximum reaction time to be simulated (see details). Default: 15.
#'
#' @return Returns a `data.frame` with three columns and `n` rows. Column names are `rt` (response
#' time), `response` (1 or 2, indicating which accumulator hit its boundary first), and `conf` (the
#' value of the confidence measure; not discretized!).
#'
#'
#' @details The simulation is done by simulating discretized steps until one process reaches
#' the boundary with an update rule:
#' \deqn{\delta X_i(t) = \max (0, X_i(t) + \delta_t ((k-1)X_i(t)-\beta X_{j=i} (t) + \mu_i + \varepsilon_i (t)),}
#' with \eqn{\varepsilon_i(t) \sim N(0, (\pi \mu_i)^2 + \sigma^2 )}. If no boundary is met within the maximum time, response is
#' set to 0. After the decision, the accumulation continues for a time period (tau), until
#' the final state is used for the computation of confidence.
#'
#'
#' @author Sebastian Hellmann.
#'
#' @name rLCA
#' @importFrom stats runif
# @importFrom pracma integral
#' @aliases simulateLCA
#'
#' @examples
#' # minimal arguments
#' simus<- rLCA(n=20, mu1=1, mu2=-0.5, th1=1, th2=0.8)
#' head(simus)
#'
#' # specifying all relevant parameters
#' simus <- rLCA(n=1000, mu1 = 2.5, mu2=1, th1=1.5, th2=1.6,
#'                k=0.1, beta=0.1, SPV=0.2, tau=0.1,
#'                wx=0.8, wrt=0.2, wint=0, t0=0.2, st0=0.1,
#'                pi=0.2, sig=1)
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   if (requireNamespace("MASS", quietly = TRUE)) {
#'     require(MASS)
#'     require(ggplot2)
#'     ggplot(simus, aes(x=rt, y=conf))+
#'       geom_bin2d()+
#'       facet_wrap(~response)
#'   }
#' }
#' boxplot(conf~response, data=simus)
#'

## When given vectorised parameters, n is the number of replicates for each parameter set
#' @rdname rLCA
#' @export
rLCA <- function (n, mu1, mu2, th1, th2,
                  k=0, beta=0, SPV=0,tau=0,
                  wx=1, wrt=0, wint=0, t0=0, st0=0,
                  pi=0, sig=1, time_scaled=TRUE, simult_conf=FALSE,
                  delta=0.01, maxrt=15)
{

  if (!time_scaled) {
    wint <- 0
    wrt <- 0
    wx <- 1
  }
  pars <- prepare_LCA_parameter(n, mu1, mu2, th1, th2, k, beta, SPV,tau,
                                wx, wrt, wint, t0, st0, pi, sig)
  res <- matrix(NA, nrow=n, ncol=6)
  for (i in seq_len(length(pars$parameter_indices))) {
    ok_rows <- pars$parameter_indices[[i]]
    current_n <- length(ok_rows)
    out <- r_LCA(current_n, pars$params[ok_rows[1], 1:15],
                 delta = delta, maxT =maxrt)

    ws <- pars$params[ok_rows[1],9:11]
    ws <- ws/sum(ws)

    # Since the Cpp function uses a different parametrization, -out[,3] is
    # exactly the distance of the loosing accumulator from its boundary
    if (time_scaled) {
      res[ok_rows,6] <- -ws[1]*out[,3] + ws[2]/sqrt(out[,1]) - ws[3]*out[,3]/sqrt(out[,1])
    } else {
      res[ok_rows,6] <- -out[,3]
    }
    # Add non-decision time to response time
    out[,1] <- out[,1] + pars$params[ok_rows[1],12] + runif(current_n, 0, pars$params[ok_rows[1],13])

    # Add inter-rating time if confidence was reported simultaneously with decision
    if (simult_conf) {
      out[1,] <- out[1,] + pars$params[ok_rows[1], 8]
    }
    res[ok_rows,1:5] <- out
  }
  res <- as.data.frame(res)
  names(res) <- c("rt", "response", "xl", "x1", "x2", "conf")
  return(res)
}




prepare_LCA_parameter <- function(nn, mu1, mu2, th1, th2,
                                  k, beta, SPV,tau,
                                  wx, wrt, wint, t0, st0,
                                  pi, sig) {
  if(any(missing(mu1), missing(mu2), missing(th1), missing(th2))) stop("mu1, mu2, th1, and th2 must be supplied")
  if ( (length(mu1) == 1) &
       (length(mu2) == 1) &
       (length(th1) == 1) &
       (length(th2) == 1) &
       (length(k) == 1) &
       (length(beta) == 1) &
       (length(SPV) == 1) &
       (length(tau) == 1) &
       (length(wx) == 1)&
       (length(wrt) ==1) &
       (length(wint) ==1) &
       (length(t0) == 1) &
       (length(st0) ==1) &
       (length(pi) ==1) &
       (length(sig) ==1)) {
    skip_checks <- TRUE
  } else {
    skip_checks <- FALSE
  }

  if (!skip_checks) {
    # all parameters brought to length of n
    mu1 <- rep(mu1, length.out = nn)
    mu2 <- rep(mu2, length.out = nn)
    th1 <- rep(th1, length.out = nn)
    th2 <- rep(th2, length.out = nn)
    k <- rep(k, length.out = nn)
    beta <- rep(beta, length.out = nn)
    SPV <- rep(SPV, length.out = nn)
    tau <- rep(tau, length.out = nn)
    wx <- rep(wx, length.out = nn)
    wrt <- rep(wrt, length.out = nn)
    wint <- rep(wint, length.out = nn)
    t0 <- rep(t0, length.out = nn)
    st0 <- rep(st0, length.out = nn)
    pi <- rep(pi, length.out = nn)
    sig <- rep(sig, length.out = nn)
  }

  # Build parameter matrix (and divide a, v, sv, sigvis and by s )

  params <- cbind (mu1, mu2, th1, th2, k, beta, SPV,tau,
                    wx, wrt, wint, t0, st0, pi, sig)

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
    parameter_indices <- list(
      seq_len(nn)
    )
  }
  list(
    params = params
    , parameter_indices = parameter_indices
  )
}

