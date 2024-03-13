#' Simulation of confidence ratings and RTs in dynaViTE with various manipulations
#'
#' Simulates the decision responses and reaction times together with a
#' discrete confidence judgment in the dynaViTE model (the 2DSD model (Pleskac & Busemeyer, 2010)
#' and the dynWEV model (Hellmann et al., 2023) are special cases), given
#' previously fitted parameters together with manipulations.
#' For details about the parameters, see \code{\link{dWEV}}.
#'
#' @param beta named vector with fitted parameters. Parameters manipulated experimentally
#' should have a parameter for each predictor named in the form "*parameter*_*predictor*".
#' The first argument can also be provided as a list with elements "beta", "fixed",
#' "model_matrix", "maxt0", and "restr_tau" providing the values for the respective arguments.
#' @param fixed list or named vector providing the parameters, which were not fitted but fixed in the model fitting.
#' @param pred_matrix matrix providing the values for the predictors used for simulation. If NULL (default),
#' the arguments model_matrix and n have to be specified for then n rows (with replacement) are drawn from
#' the model_matrix are drawn for prediction.
#' @param model_matrix matrix providing the values of the predictors used for model fitting (only necessary, if pred_matrix is NULL).
#' @param n integer. The number of samples drawn from the model_matrix (only necessary, if pred_matrix is NULL).
#' @param maxt0 numeric. Maximum possible value for t0 in the fitting procedure. t0 represents minimal non-decision time,
#' so it is bounded from above by the minimal observed response time in the data. This is used to scale the fitted parameter for t0 properly.
#' @param restr_tau numeric or Inf. Maximum possible value for tau in the fitting procedure. tau represents post-decisional accumulation time,
#' and my be bounded from above in the fitting (e.g. based on observed confidence response times).
#' This is used to scale the fitted parameter for t0 properly.
#' @param delta numeric. Discretization steps for simulations with the stochastic process.
#' @param maxrt numeric. Maximum reaction time returned.
#' If the simulation of the stochastic process exceeds a rt of `maxrt`,
#' the response will be set to 0 and `maxrt` will be returned as rt.
#' @param  simult_conf logical. `TRUE`, if in the experiment confidence was reported simultaneously
#' with the decision, as then decision and confidence judgment are assumed to have happened
#' subsequent before response and tau is added to the simulated decision time. If `FALSE`
#' returned response time will only be decision time plus non-judgment time component.
#'
#' @param seed numerical. Seeding for non-random data generation.
#' @param process_results logical. Whether the output simulations should contain the final
#' state of the decision (and visibility) process as additional column. Default is FALSE, meaning that
#' no additional columns for the final process states are returned.
#'
#' @return Returns a `matrix` with columns for each predictor in the model_matrix/pred_matrix and
#' simulated outcome variables: `response`, `rt`, `conf` (the continuous confidence
#' measure) and `rating` (the discrete confidence rating), and `dec` and `vis`
#' (only if `process_results=TRUE`) for the final states of accumulators in the
#' simulation.
#'
#'
#' @details Simulation of responses and decision times is done by simulating
#' normal variables in discretized steps until the lower or upper boundary
#' is met (or the maximal rt is reached). Afterwards, a confidence measure
#' is simulated according to the respective model.
#'
#' The confidence outputs are then binned according to the given thresholds.
#' The output of the fitting function \code{\link{fitRTConf_formula}} fits the argument `beta` for simulation.
#'
#' @note Different parameters for different conditions are only allowed for drift rate,
#' \code{v}, drift rate variability, \code{sv} and diffusion constant `s`.
#' All other parameters are used for all conditions.
#'
#' @references  Hellmann, S., Zehetleitner, M., & Rausch, M. (2023). Simultaneous modeling of choice, confidence and response time in visual perception. \emph{Psychological Review} 2023 Mar 13. doi: 10.1037/rev0000411. Epub ahead of print. PMID: 36913292.
#'
#'
#' @author Sebastian Hellmann.
#'
#' @name simulateWEV_formula
#'
#' @examples
#' # Examples for "dynWEV" model (equivalent applicable
#' # for "2DSD" model (with less parameters))
#' # 1. Define some parameter set in a data.frame
#' paramDf <- data.frame(a=2.5,v1=0.1, v2=1, t0=0.1,z=0.55,
#'                       sz=0.3,sv=0.8, st0=0,  tau=3, w=0.1,
#'                       theta1=0.8, svis=0.5, sigvis=0.8)
#'
#' # 2. Simulate trials for both stimulus categories and all conditions (2)
#' simus <- simulateWEV(paramDf, model="dynWEV")
#' head(simus)
#' \donttest{
#'   library(ggplot2)
#'   simus <- simus[simus$response!=0,]
#'   simus$rating <- factor(simus$rating, labels=c("unsure", "sure"))
#'   ggplot(simus, aes(x=rt, group=interaction(correct, rating),
#'                     color=as.factor(correct), linetype=rating))+
#'     geom_density(linewidth=1.2)+xlim(c(0,5))+
#'     facet_grid(rows=vars(condition), labeller = "label_both")
#' }
#'
#' # automatically aggregate simulation distribution
#' # to get only accuracy x confidence rating distribution for
#' # all conditions
#' agg_simus <- simulateWEV(paramDf, model="dynWEV", agg_simus = TRUE)
#' head(agg_simus)
#' \donttest{
#'   agg_simus$rating <- factor(agg_simus$rating, labels=c("unsure", "sure"))
#'   library(ggplot2)
#'   ggplot(agg_simus, aes(x=rating, group=correct, fill=as.factor(correct), y=p))+
#'     geom_bar(stat="identity", position="dodge")+
#'     facet_grid(cols=vars(condition), labeller = "label_both")
#' }
#' \donttest{
#'   # Compute Gamma correlation coefficients between
#'   # confidence and other behavioral measures
#'   # output will be a list
#'   simu_list <- simulateWEV(paramDf,n = 400, model="dynWEV", gamma=TRUE)
#'   simu_list
#' }

#' @rdname simulateWEV_formula
#' @export
simulate_dynaViTE_formula <- function(beta, pred_matrix=NULL, fixed=NULL, #method="simulation_matrix",
                                     model_matrix=NULL, n=NULL, maxt0=NULL, restr_tau=NULL,
                                     delta=0.01, maxrt=15, simult_conf = FALSE,
                                     process_results=FALSE, seed=NULL) {
  if (is.list(beta) && all(c("beta", "model_matrix", "fixed", "restr_tau", "maxt0") %in% names(beta))) {
    model_matrix <- beta$model_matrix
    fixed <- beta$fixed
    maxt0 <- beta$maxt0
    restr_tau <- beta$restr_tau
    beta <- beta$beta
  }

  # turn 'fixed' from character to a list
  if (is.character(fixed)) eval(parse(text=paste0("fixed <- list(",fixed,")")))
  #betatrans <- transform_const_pars_to_scale(beta, fixed, maxt0=res_22$maxt0)

  manipulations <- list()
  man_pars <- grep(pattern="_", names(beta), value=TRUE)
  man_pars <- strsplit(man_pars, "_")
  for (i in 1:length(man_pars)) {
    manipulations[[man_pars[[i]][1]]] <- c(manipulations[[man_pars[[i]][1]]],
                                           man_pars[[i]][2])
  }


  # if (method=="simulation_matrix") {     ## So far the only method
  if (is.null(pred_matrix)) {
    if (is.null(n)) stop("Please specify the n argument")
    pred_matrix <- model_matrix[sample(1:nrow(model_matrix), replace=TRUE, size = n), ]
  }

  #parnames <- c("v", "z", "a", "sz", "t0", "st0", "d", "sv", "tau", "w", "svis", "sigvis", "lambda", "s", "th1", "th2")
  parnames <- c("a", "v", "t0", "d","z",  "sz", "sv", "st0","tau", "th1", "th2", "lambda", "w", "muvis", "sigvis", "svis", "s")
  fixed01 <- c("z", "sz", "w", "d", "st0")
  fixedpos <- c("v", "a", "t0", "sv", "svis", "sigvis", "lambda", "s")
  parammatrix <-  matrix(NA, nrow=nrow(pred_matrix), ncol=length(parnames))
  colnames(parammatrix) <- parnames
  for (i in 1:length(parnames)) {
    if (parnames[i] %in% names(fixed)) {   ### If parameter was fixed, use fixed value for column (no transformation, then!)
      parammatrix[,parnames[i]] <- fixed[[parnames[i]]]
    } else {
      # If parameter was not fixed, then either it was manipulated:
      if (parnames[i] %in% names(manipulations)) {
        parammatrix[,parnames[i]] <-
          pred_matrix[,manipulations[[parnames[i]]], drop=FALSE] %*% # take columns from prediction model matrix and multiply
          beta[grepl(names(beta), pattern=paste0(parnames[i], "_"))] # with respective beta-values
      } else { # or it was fit as a constant, then take parameter as fitted:
        parammatrix[,parnames[i]] <- beta[parnames[i]]
      }
      ## Afterwards, for all parameters involved in the fitting procedure,
      # we have to transform them back to their model scale:
      if (parnames[i] %in% fixed01) parammatrix[,parnames[i]] <- pnorm(parammatrix[,parnames[i]])
      if (parnames[i] %in% fixedpos) parammatrix[,parnames[i]] <- exp(parammatrix[,parnames[i]])
      if (parnames[i] == "sz") parammatrix[,"sz"] <- (pmin(parammatrix[,"z"], (1-parammatrix[,"z"]))*2)*parammatrix[,"sz"]
      if (parnames[i] == "t0") parammatrix[,"t0"] <- maxt0*parammatrix[,"t0"]
      if (parnames[i] == "d") parammatrix[,"d"] <- parammatrix[,"t0"]*parammatrix[,"d"]
      if (parnames[i] == "tau") {
        if (restr_tau == Inf) {
          parammatrix[,parnames[i]] <- exp(parammatrix[,parnames[i]])
        } else if (simult_conf) {
          parammatrix[,parnames[i]] <- pnorm(parammatrix[,parnames[i]])*(maxt0-parammatrix[,"t0"])
        } else {
          parammatrix[,parnames[i]] <- restr_tau * pnorm(parammatrix[,parnames[i]])
        }
      }
    }

  }
  # Fill single missing muvis (index: 14) values with absolute value of decision drift v
  parammatrix[is.na(parammatrix[,14]), 14] <-
    abs(parammatrix[is.na(parammatrix[,14]), 2])

  # Scale the parameters scaled by s by s
  # In the actual function:
  # cbind (a/s, v/s, t0, d, sz, sv/s, st0, z,
  #        tau, th1/s, th2/s, lambda, w, muvis/s, sigvis/s, svis/s, numeric_bounds)
  parammatrix[, c(1, 2, 7, 10, 11, 14, 15, 16)] <-
    parammatrix[, c(1, 2,7, 10, 11, 14, 15, 16)]/parammatrix[,17]

  # bound between trial variability st0 to 2 (seconds)
  parammatrix[,8] <-  parammatrix[,8]*2


  # t0 <- t0+st0/2
  parammatrix[,3] <- parammatrix[,3] + parammatrix[,8]/2

  ## Simulation:
  if (!is.null(seed)) set.seed(seed)
  res <- r_WEV_matrix(parammatrix[,c(1:4, 6:8, 5, 9, 12:16)], delta=delta,
                      maxT = maxrt, stop_on_error = TRUE)
  if (process_results) {
    colnames(res) <- c("rt", "response", "conf", "dec", "vis", "mu")
  } else {
    res <- res[,c(1:3)]
    colnames(res) <- c("rt", "response", "conf")
  }

  # Compute confidence rating:
  thetas <- beta[grep(x=names(beta), pattern="theta")]
  thetasUpper <- thetas[grep(x=names(thetas), pattern="Upper")]
  thetasLower <- thetas[grep(x=names(thetas), pattern="Lower")]
  if (length(thetasLower)==0 && length(thetasUpper)==0) {
    thetasLower <- thetas
    thetasUpper <- thetas
  }
  nRatings <- length(thetasUpper)+1
  thetasUpper <- c(-Inf,thetasUpper, Inf)
  thetasLower <- c(-Inf,thetasLower, Inf)


  levels_lower <- cumsum(as.numeric(table(thetasLower)))
  levels_lower <- levels_lower[-length(levels_lower)]
  levels_upper <- cumsum(as.numeric(table(thetasUpper)))
  levels_upper <- levels_upper[-length(levels_upper)]
  thetasLower <- unique(thetasLower)
  thetasUpper <- unique(thetasUpper)

  res <- cbind(res,1)
  res[res[,2]==1, 4] <- as.numeric(as.character(cut(res[res[,2]==1, 3],
                                                    breaks=thetasUpper, labels = levels_upper)))
  res[res[,2]==-1,4] <- as.numeric(as.character(cut(res[res[,2]==-1, 3],
                                                    breaks=thetasLower, labels = levels_lower)))

  if (simult_conf) {
    res[,1] <- res[,1] + parammatrix[,'tau']
  }
  colnames(res)[4] <- "rating"
  res <- cbind(pred_matrix, res)
  return(res)
}

#' @rdname simulateWEV_formula
#' @export
transform_const_pars_to_scale <- function(beta, fixed=list(), maxt0=NULL, restr_tau=Inf, simult_conf=FALSE) {
  parnames <- c("a", "v", "t0", "d","z",  "sz", "sv", "st0","tau", "th1", "th2", "lambda", "w", "muvis", "sigvis", "svis", "s")
  fixed01 <- c("z", "sz", "w", "d", "st0")
  fixedpos <- c("v", "a", "t0", "sv", "svis", "sigvis", "lambda", "s")
  for (i in 1:length(parnames)) {
    if (parnames[i] %in% names(fixed)) {
      beta[parnames[i]] <- fixed[[parnames[i]]]
    } else if (parnames[i] %in% names(beta)) {
      if (parnames[i] %in% fixed01) beta[parnames[i]] <- pnorm(beta[parnames[i]])
      if (parnames[i] %in% fixedpos) beta[parnames[i]] <- exp(beta[parnames[i]])
      if (parnames[i] == "sz") beta["sz"] <- (pmin(beta["z"], (1-beta["z"]))*2)*beta["sz"]
      if (parnames[i] == "t0") beta["t0"] <- maxt0*beta["t0"]
      if (parnames[i] == "d") beta["d"] <- beta["t0"]*beta["d"]
      if (parnames[i] == "tau") {
        if (restr_tau == Inf) {
          beta[parnames[i]] <- exp(beta[parnames[i]])
        } else if (simult_conf) {
          beta[parnames[i]] <- pnorm(beta[parnames[i]])*(maxt0-beta["t0"])
        } else {
          beta[parnames[i]] <- restr_tau * pnorm(beta[parnames[i]])
        }
      }
    }

  }
  return(beta)
}
