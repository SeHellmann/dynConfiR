#' Function for fitting several sequential sampling confidence models in parallel
#'
#' This function is a wrapper of the function \code{\link{fitConfModel}} (see
#' there for more information). It calls the function for every possible combination
#' of model and participant in `model` and \code{data} respectively.
#' Also, see \code{\link{dWEV}}, \code{\link{d2DSD}}, \code{\link{dDDMConf}},
#' and \code{\link{dRM}} for more
#' information about the parameters.
#'
#' @param data a `data.frame` where each row is one trial, containing following
#' variables (column names can be changed by passing additional arguments of
#' the form \code{condition="contrast"}):
#' * \code{condition} (not necessary; for different levels of stimulus quality, will be transformed to a factor),
#' * \code{rating} (discrete confidence judgments, should be given as integer vector; otherwise will be transformed to integer),
#' * \code{rt} (giving the reaction times for the decision task),
#' * either 2 of the following (see details for more information about the accepted formats):
#'   * \code{stimulus} (encoding the stimulus category in a binary choice task),
#'   * \code{response} (encoding the decision response),
#'   * \code{correct} (encoding whether the decision was correct; values in 0, 1)
#' * \code{sbj} (giving the subject ID; the models given in the second argument are fitted for each
#'   subject individually. (Furthermore, if `logging = TRUE`, the ID is used in files
#'   saved with interim results and logging messages.))
#' @param models character vector with following possible elements "dynWEV", "2DSD", "IRM", "PCRM", "IRMt", and "PCRMt"  for the models to be fit.
#' @param nRatings integer. Number of rating categories. If `NULL`, the maximum of
#' `rating` and `length(unique(rating))` is used. This argument is especially
#' important for data sets where not the whole range of rating categories is realized.
#' If given, ratings has to be given as factor or integer.
#' @param fixed list. List with parameter value pairs for parameters that should not be fitted. (see Details).
#' @param restr_tau numerical or `Inf` or `"simult_conf"`. Used for 2DSD and dynWEV only. Upper bound for tau.
#' Fits will be in the interval (0,`restr_tau`). If `FALSE` tau will be unbound. For `"simult_conf"`, see the documentation of
#' \code{\link{d2DSD}} and \code{\link{dWEV}}
#' @param grid_search logical. If `FALSE`, the grid search before the optimization
#' algorithm is omitted. The fitting is then started with a mean parameter set
#' from the default grid. (Default: `TRUE`)
#' @param opts list. A list for more control options in the optimization routines
#' (depending on the `optim_method`). See details for more information.
#' @param optim_method character. Determines which optimization function is used for
#' the parameter estimation. Either `"bobyqa"` (default), `"L-BFGS-B"` or `"Nelder-Mead"`.
#' `"bobyqa"` uses a box-constrained optimization with quadratic interpolation.
#' (See \code{\link[minqa]{bobyqa}} for more information.) The first two use a
#' box-constraint optimization. For Nelder-Mead a transfinite function rescaling is used
#' (i.e. the constrained arguments are suitably transformed to the whole real line).
#' @param logging logical. If `TRUE`, a folder 'autosave/fit**model**' is created and
#' messages about the process are printed in a logging file and to console (depending
#' on OS). Additionally intermediate results are saved in a `.RData` file with the
#' participant ID in the name.
#' @param parallel "models", "single", "both" or `FALSE`. If `FALSE` no parallelization
#' is used in the fitting process. If "models" the fitting process is parallelized over
#' participants and models (i.e. over the calls for fitting functions). If "single"
#' parallelization is used within the fitting processes (over initial grid search and
#' optimization processes for different start points, but see \code{\link{fitRTConf}}).
#' If "both", parallelization is done hierarchical. For small number of
#' models and participants "single" or "both" is preferable. Otherwise, you may use "models".
#' @param precision numerical scalar. For 2DSD and dynWEV only. Precision of calculation.
#' (in the respective models) for the density functions (see \code{\link{dWEV}} for more information).
#' @param n.cores integer vector or `NULL`. If \code{parallel} is "models" or "single", a single
#' integer for the number of cores used for parallelization is required. If
#' \code{parallel} is "both", two values are required. The first for the number of parallel
#' model-participant combinations and the second for the parallel processes within the
#' fitting procedures (this may be specified
#' to match the \code{nAttemps}-Value in the \code{opts} argument. If `NULL` (default)
#' the number of available cores -1 is used.
#' If `NULL` and \code{parallel} is "both", the cores will be used for
#' model-participant-parallelization, only.
#' @param ... Possibility of giving alternative variable names in data frame
#' (in the form \code{condition = "SOA"}, or \code{response="pressedKey"}).
#'
#' @return Gives data frame with rows for each model-participant combination and columns for the different parameters
#' as fitted result as well as additional information about the fit (`negLogLik` (for final parameters),
#' `k` (number of parameters), `N` (number of data rows), `BIC`, `AICc` and `AIC`)
#'
#' @details The fitting involves a first grid search through an initial grid. Then the best \code{nAttempts}
#' parameter sets are chosen for an optimization, which is done with an algorithm, depending on the argument
#' \code{optim-method}. The Nelder-Mead algorithm uses the R function \code{\link[stats]{optim}}.
#' The optimization routine is restarted \code{nRestarts} times with the starting parameter set equal to the
#' best parameters from the previous routine.
#'
#'  \strong{stimulus, response and correct}. Two of these columns must be given in data. If all three are given, correct will have no effect (and will be not checked!).
#'  stimulus can always be given in numerical format with values -1 and 1. response can always be given as a character vector with "lower" and "upper" as values.
#'  Correct must always be given as a 0-1-vector. If stimulus is given together with response and they both do not match the above format, they need to have the same values/levels (if factor).
#'  In the case that only stimulus/response is given in any other format together with correct, the unique values will be sorted increasingly and
#'  the first value will be encoded as "lower"/-1 and the second as "upper"/+1.
#'
#'  \strong{fixed}. Parameters that should not be fitted but kept constant. These will be dropped from the initial grid search
#'  but will be present in the output, to keep all parameters for prediction in the result. Includes the
#' possibility for symmetric confidence thresholds for both alternative (\code{sym_thetas}=logical). Other examples are
#' \code{z =.5}, \code{sv=0}, \code{st0=0}, \code{sz=0}. For race models, the possibility of setting \code{a='b'} (or vice versa)
#' leads to identical upper bounds on the decision processes, which is the equivalence for \code{z=.5} in a diffusion process
#'
#'  \strong{opts}. A list with numerical values. Possible options are listed below (together with the optimization method they are used for).
#'  * \code{nAttempts} (all) number of best performing initial parameter sets used for optimization; default 5
#'  * \code{nRestarts} (all) number of successive `optim` routines for each of the starting parameter sets; default 5,
#'  * \code{maxfun} (\code{'bobyqa'}) maximum number of function evaluations; default: 5000,
#'  * \code{maxit} (\code{'Nelder-Mead' and 'L-BFGS-B'}) maximum iterations; default: 2000,
#'  * \code{reltol} (\code{'Nelder-Mead'}) relative tolerance; default:  1e-6),
#'  * \code{factr} (\code{'L-BFGS-B'}) tolerance in terms of reduction factor of the objective, default: 1e-10)
#'
#' @md
#'
#' @references  Hellmann, S., Zehetleitner, M., & Rausch, M. (2023). Simultaneous modeling of choice, confidence and response time in visual perception. \emph{Psychological Review} 2023 Mar 13. doi: 10.1037/rev0000411. Epub ahead of print. PMID: 36913292.
#'
#' @author Sebastian Hellmann, \email{sebastian.hellmann@@ku.de}
#'
#' @name fitRTConfModels_formula
#' @importFrom dplyr rename
#' @import parallel
# @importFrom pracma integral
#'
#' @examples
#' # 1. Generate data from two artificial participants
#' # Get random drift direction (i.e. stimulus category) and
#' # stimulus discriminability (two steps: hard, easy)
#' stimulus <- sample(c(-1, 1), 400, replace=TRUE)
#' discriminability <- sample(c(1, 2), 400, replace=TRUE)
#'
#' # generate data for participant 1
#' data <- rWEV(400, a=2, v=stimulus*discriminability*0.5,
#'              t0=0.2, z=0.5, sz=0.1, sv=0.1, st0=0,  tau=4, s=1, w=0.3)
#' # discretize confidence ratings (only 2 steps: unsure vs. sure)
#' data$rating <- as.numeric(cut(data$conf, breaks = c(-Inf, 1, Inf), include.lowest = TRUE))
#' data$participant = 1
#' data$stimulus <- stimulus
#' data$discriminability <- discriminability
#' # generate data for participant 2
#' data2 <- rWEV(400, a=2.5, v=stimulus*discriminability*0.7,
#'              t0=0.1, z=0.7, sz=0, sv=0.2, st0=0,  tau=2, s=1, w=0.5)
#' data2$rating <- as.numeric(cut(data$conf, breaks = c(-Inf, 0.3, Inf), include.lowest = TRUE))
#' data2$participant = 2
#' data2$stimulus <- stimulus
#' data2$discriminability <- discriminability
#'
#' # bind data from participants
#' data <- rbind(data, data2)
#' data <- data[data$response!=0, ] # drop not finished decision processes
#' data <- data[,-3] # drop conf measure (unobservable variable)
#' head(data)
#'
#'
#' # 2. Use fitting function
#' \dontrun{
#'   # Fitting takes very long to run and uses multiple (6) cores with this
#'   # call:
#'   fitRTConfModels(data, models=c("dynWEV", "PCRM"), nRatings = 2,
#'                 logging=FALSE, parallel="both",
#'                 n.cores = c(2,3), # fit two participant-model combination in parallel
#'                 condition="discriminability")# tell which column is "condition"
#' }
#'
#'


#' @rdname fitRTConfModels_formula
#' @export
fitRTConfModels_formula <- function(data, models = c("dynaViTE",
                                                     "PCRMt", "IRMt"),
                                    manipulations=NULL,
                                     fixed = list(sym_thetas = FALSE), nRatings = NULL,
                                    restr_tau=Inf, grid_search=TRUE,
                                    opts=list(), optim_method = "Nelder-Mead", logging=FALSE, precision=2,
                                    parallel = TRUE, n.cores=NULL, ...){ #  ?ToDO: vary_sv=FALSE, RRT=NULL, vary_tau=FALSE
  if (any(!grepl("dynaViTE|IRM|PCRM|IRMt|PCRMt|dynWEV|2DSD|DDMConf", models))) {
    stop("all models must be one of:
    'dynaViTE', 'dynWEV', '2DSDT', '2DSD',
    'IRM', 'IRMt', or 'PCRMt'")
  }
  if (any(grepl("IRM|PCRM|IRMt|PCRMt|DDMConf", models))) stop("Only dynaViTE (and 2DSD, and dynWEV) implemented so far!")
  ### Maybe later: use ...-argument für renaming data-columns and to pass other arguments
  # colrenames <- c(...)
  # if (length(colrenames)>0) {
  #   addargs <- list(colrenames[!(colrenames %in% names(data))])
  #   list2env(addargs, envir=environment())
  #   colrenames <- colrenames[colrenames %in% names(data)]
  #   data <- rename(data, colrenames)
  # }
  tryCatch(data <- rename(data, ...),
           error = function(e) stop("Error renaming data columns. Probably a column name does not exist.\nCheck whether an argument was misspelled and data name pairs are given in the form expected_name = true_name."))

  ### Adapt arguments for individual settings
  if (!("sym_thetas" %in% names(fixed))) fixed["sym_thetas"] <- FALSE
  sym_thetas <- as.logical(fixed['sym_thetas'])


  ### Determine number of jobs, i.e. model-participant-combinations
  sbjcol <- c("subject", "participant", "sbj")[which(c("subject", "participant", "sbj") %in% names(data))]
  if (length(sbjcol)==0) {
    data$sbj <- 999
    sbjcol <- "sbj"
  } else {
    if (sbjcol != "sbj") {
      data$sbj <- data[[sbjcol]]
      data[[sbjcol]] <- NULL
    }
  }
  subjects <- unique(data$sbj)
  nJobs <- length(models)*length(subjects)

  ### set up parallelization settings
  parallel.models <- FALSE
  parallel.single <- FALSE
  n.cores.models = 0
  n.cores.single = 0
  if (parallel == "both") {
    parallel.models <- TRUE
    parallel.single <- TRUE
    if (is.null(n.cores)) {
      n.cores.models <- min(detectCores()-1, nJobs)
      parallel.single <- FALSE
    } else {
      n.cores.models <- n.cores[1]
      n.cores.single <- n.cores[2]
    }
  } else if (parallel == "single") {
    parallel.single <- TRUE
    if (is.null(n.cores)) {
      n.cores.single <- detectCores()-1
    } else {
      n.cores.single <- n.cores
    }
  } else if (parallel == "models") {
    parallel.models <- TRUE
    if (is.null(n.cores)) {
      n.cores.models <- min(detectCores()-1, nJobs)
    } else {
      n.cores.models <- n.cores
    }
  }



  nConds <- length(unique(data$condition))
  rating <- data$rating
  if (!is.numeric(rating)){
    rating <- as.integer(as.factor(rating))
  }
  if (!is.integer(rating)){
    rating <- as.integer(rating)
  }
  if (is.integer(rating)) {
    if (is.null(nRatings)) {
      nRatings <- max(rating)
    }
    if (max(length(unique(rating)), max(rating)+1-min(rating))>nRatings && any(rating==0)) {
      rating <- rating + 1L
      nRatings <- nRatings + 1
    }
  }
  if (!inherits(manipulations[[1]], "formula") &
      ((length(models) != length(manipulations)) ||
       !inherits(manipulations[[1]][[1]], "formula"))) {
    stop("models and manipulations arguments do not match.
         Either provide several models and  one list of manipulations with formulas as elements or
         a list of manipulation-lists with the same length as the models argument.")
  }

  call_fitfct <- function(X) {
    if (inherits(manipulations[[1]], "formula")) {
      cur_manipulations <- manipulations
    } else {
      cur_manipulations <- manipulations[[X[[1]]]]
    }
    cur_model <- models[X[[1]]]
    cur_sbj <- X[[2]]
    sbj <- NULL # to omit a note because of an unbound variable
    data_part <- subset(data, sbj==cur_sbj)

    res <- fitRTConf_formula(data_part, model = cur_model,
                             manipulations=cur_manipulations,
                             fixed = fixed, nRatings = nRatings,
                             restr_tau =restr_tau,
                             grid_search = grid_search, precision=precision,
                             logging=logging, opts=opts, optim_method = optim_method,
                             useparallel = parallel.single, n.cores=n.cores.single)
    res$model <- cur_model
    res$sbj <- cur_sbj
    # res[outnames[!(outnames %in% names(res))]] <- NA
    # res <- res[,outnames]
    return(res)
  }

  jobs <- expand.grid(model=1:length(models), sbj=subjects)
  listjobs <- list()
  for (i in 1:nrow(jobs)) {
    listjobs[[i]] <- list(model = jobs[["model"]][i], sbj = jobs[["sbj"]][i])
  }
  if (parallel.models) {

    clmodels <- makeCluster(type="SOCK", n.cores.models)
    clusterExport(clmodels, c("data", "manipulations", "fixed", "nRatings", "restr_tau", "models",
                              "logging", "opts", "optim_method", "grid_search",
                              "parallel.single", "n.cores.single", "precision",
                              "call_fitfct"), envir = environment())
    on.exit(try(stopCluster(clmodels), silent = TRUE))
    res <- clusterApplyLB(clmodels, listjobs, fun=call_fitfct)
    stopCluster(clmodels)
  } else {
    res <- lapply(X=listjobs, FUN=call_fitfct)

  }

  return(res)
}
