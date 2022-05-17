#' Function for fitting sequential sampling confidence models
#'
#' Fits the parameters of different models of response time and confidence, including
#' the 2DSD model (Pleskac & Busemeyer, 2010), dynWEV, and various
#' flavors of race models (Hellmann et al., preprint). Which model to fit is
#' specified by the argument \code{model}.
#' So far, only a ML method is implemented (a quantile-Chi square method may be
#' implemented in the future, too).
#' See \code{\link{dWEV}}, \code{\link{d2DSD}}, and \code{\link{dRM}} for more
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
#' * \code{sbj} or \code{participant} (optional; giving the subject ID; only relevant if logging == TRUE;
#'                                                       if unique the ID is used in autosave files and logging messages;
#'                                                       if non-unique or missing AND logging ==TRUE, 999 will be used then)
#' @param model character scalar. One of "dynWEV", "2DSD", "IRM", "PCRM", "IRMt", "PCRMt" for the model to be fit.
#' @param fixed list. List with parameter value pairs for parameters that should not be fitted. See Details.
#' @param init_grid data.frame or `NULL`. Grid for the initial parameter search. Each row is one parameter constellation.
#' See details for more information. If \code{NULL} a default grid will be used.
#' @param grid_search logical. If `FALSE`, the grid search before the optimization
#' algorithm is omitted. The fitting is then started with a mean parameter set
#' from the default grid (if `init_grid=NULL`) or directly with the rows from
#' `init_grid`, if not NULL. (Default: `TRUE`)
#' @param data_names named list (e.g. `c(rating="confidence")`). Alternative
#' possibility of giving other column names for the variables in the data. By default
#' column names are identical to the ones given in the data argument description.
#' @param nRatings integer. Number of rating categories. If `NULL`, the maximum of
#' `rating` and `length(unique(rating))` is used. This argument is especially
#' important for data sets where not the whole range of rating categories is realized.
#' If given, ratings has to be given as factor or integer.
#' @param logging logical. If `TRUE`, a folder autosave/fit**model** is created and
#' messages about the process are printed in a logging file and to console (depending
#' on OS). Additionally intermediate results are saved in a .RData file with the
#' participant ID in the name.
#' @param opts list. A list for more control options in the optimization routines
#' (depending on the `optim_method`). See details for more information.
#' @param optim_method character. Determines which optimization function is used for
#' the parameter estimation. Either `"bobyqa"` (default), `"L-BFGS-B"` or `"Nelder-Mead"`.
#' `"bobyqa"` uses a box-constrained optimization with quadratic interpolation.
#' (See \code{\link[minqa]{bobyqa}} for more information.) The first two use a
#' box-constraint optimization. For Nelder-Mead a transfinite function rescaling is used
#' (i.e. the constrained arguments are suitably transformed to the whole real line).
#' @param useparallel logical. If `TRUE` the grid search in the beginning is done with a
#' parallel back-end, using the \code{parallel} package.
#' @param n.cores integer or `NULL`. Number of cores used for parallelization. If `NULL`
#' (default) the number of available cores -1 is used.
#' @param restr_tau numerical or Inf or "simult_conf". For 2DSD and dynWEV only.
#' Upper bound for tau. Fits will be in the interval (0,restr_tau). If FALSE tau will be unbound.
#' For "simult_conf", see the documentation of \code{\link{d2DSD}} and \code{\link{dWEV}}
#' @param precision numerical scalar. For 2DSD and dynWEV only. Precision of calculation.
#' (in the respective models) for the density functions (see \code{\link{dWEV}} for more information).
#' @param ... Possibility of giving alternative variable names in data frame
#' (in the form \code{condition = "SOA"}, or \code{response="pressedKey"}).
#'
#' @return Gives a one-row data frame with columns for the different parameters as
#' fitted result as well as additional information about the fit (`negLogLik` (for
#' final parameters), `k` (number of parameters), `N` (number of data rows),
#' `BIC`, `AICc` and `AIC`) and the column `fixed`, which includes all information
#' about fixed and not fitted parameters.
#'
#' @details The fitting involves a first grid search through computation of the
#' likelihood on an initial grid with possible sets of parameters to start the
#' optimization routine. Then the best \code{nAttempts} parameter sets are
#' chosen for an optimization, which is done with an algorithm, depending on the
#' argument \code{optim-method}. The Nelder-Mead algorithm uses the R function
#' \code{\link[stats]{optim}}. The optimization routine is restarted
#' \code{nRestarts} times with the starting parameter set equal to the
#' best parameters from the previous routine.
#'
#'  \strong{stimulus, response and correct}. Two of these columns must be given in
#'  data. If all three are given, correct will have no effect (and will be not checked!).
#'  stimulus can always be given in numerical format with values -1 and 1. response
#'  can always be given as a character vector with "lower" and "upper" as values.
#'  Correct must always be given as a 0-1-vector. If stimulus is given together with
#'  response and they both do not match the above format, they need to have the same
#'  values/levels (if factor).
#'  In the case that only stimulus/response is given in any other format together with
#'  correct, the unique values will be sorted increasingly and
#'  the first value will be encoded as "lower"/-1 and the second as "upper"/+1.
#'
#'  \strong{fixed}. Parameters that should not be fitted but kept constant. These will
#'  be dropped from the initial grid search
#'  but will be present in the output, to keep all parameters for prediction in the result.
#'  Includes the possibility for symmetric confidence thresholds for both alternative
#'  (\code{sym_thetas}=logical). Other examples are
#' \code{z =.5}, \code{sv=0}, \code{st0=0}, \code{sz=0}. For race models, the possibility
#' of setting \code{a='b'} (or vice versa)
#' leads to identical upper bounds on the decision processes, which is the equivalence for
#'  \code{z=.5} in a diffusion process
#'
#' \strong{init_grid}. Each row should be one parameter set to check. The column names
#' should include the parameters of the desired model, which are the following for 2DSD:
#' a, vmin and vmax (will be equidistantly spanned across conditions), sv, z (as the
#' relative starting point between 0 and a), sz (also in relative terms), t0, st0, theta0
#' (minimal threshold), thetamax (maximal threshold; the others will be equidistantly
#' spanned symmetrically for both decisions), and tau. For dynWEV,
#' additionally w , svis, and sigvis are required. For the race models the parameters
#' are: vmin, vmax (will be equidistantly
#' spanned across conditions), a and b (decision thresholds), t0, st0, theta0 (minimal
#'  threshold), thetamax (maximal threshold;
#' the others will be equidistantly spanned symmetrically for both decisions), and for
#' time-dependent confidence race models
#' additionally wrt and wint (as weights compared to wx=1).
#'
#'  \strong{opts}. A list with numerical values. Possible options are listed below
#'  (together with the optimization method they are used for).
#'  * \code{nAttempts} (all) number of best performing initial parameter sets used for
#'   optimization; default 5, if grid_search is TRUE.
#'  If grid_search is FALSE and init_grid is NULL, then nAttempts will be set to 1 (and
#'  any input will be ignored).
#'  If grid_search is FALSE and init_grid is not NULL, the rows of init_grid will be used
#'  from top to bottom
#'  (since no initial grid search is done) with not more than nAttempts rows used.
#'  * \code{nRestarts} (all) number of successive optim routines for each of the starting parameter sets; default 5,
#'  * \code{maxfun} (\code{'bobyqa'}) maximum number of function evaluations; default: 5000,
#'  * \code{maxit} (\code{'Nelder-Mead' and 'L-BFGS-B'}) maximum iterations; default: 2000,
#'  * \code{reltol} (\code{'Nelder-Mead'}) relative tolerance; default:  1e-6),
#'  * \code{factr} (\code{'L-BFGS-B'}) tolerance in terms of reduction factor of the objective, default: 1e-10)
#'
#' @md
#'
#' @references Hellmann, S., Zehetleitner, M., & Rausch, M. (preprint). Simultaneous modeling of choice,
#' confidence and response time in visual perception. https://osf.io/9jfqr/
#'
#' https://nashjc.wordpress.com/2016/11/10/why-optim-is-out-of-date/
#'
#' https://www.damtp.cam.ac.uk/user/na/NA_papers/NA2009_06.pdf
#'
#'
#'
#' @author Sebastian Hellmann.
#'
#' @name fitRTConf
#' @importFrom stats setNames aggregate optim qnorm pnorm
#' @importFrom minqa bobyqa
#' @importFrom dplyr if_else case_when rename
#' @import logger
#' @import parallel
#' @importFrom magrittr %>%
#' @importFrom rlang .data
# @importFrom pracma integral
#' @aliases fitSeqSampConf fitConfModel fitConf fitConfRT
#' @importFrom Rcpp evalCpp
#'

#' @rdname fitRTConf
#' @export
fitRTConf <- function(data, model = "dynWEV",
                      fixed = list(sym_thetas = FALSE),
                      init_grid = NULL, grid_search = TRUE,
                      data_names = list(), nRatings = NULL, restr_tau =Inf,
                      precision=1e-5,logging=FALSE, opts=list(), optim_method = "bobyqa",
                      useparallel = FALSE, n.cores=NULL, ...){ #  ?ToDO: vary_sv=FALSE, RRT=NULL, vary_tau=FALSE
  ### Check model argument
  if (model == "WEVmu") {  ## Old name for dynWEV model
    model <- "dynWEV"
  }
  if (!model %in% c("IRM", "PCRM", "IRMt", "PCRMt", "dynWEV", "2DSD")) {
    stop("model must be 'dynWEV', '2DSD', 'IRM', 'PCRM', 'IRMt', or 'PCRMt'")
  }

  #### Check argument types ###
  if (!is.logical(grid_search)) stop(paste("grid_search must be logical, but is ", typeof(grid_search), sep=""))
  if (length(grid_search)!=1) stop(paste("grid_search must be of length 1, it's ", length(grid_search), sep=""))

  # colrenames <- c(...)
  # colrenames <- colrenames[colrenames %in% names(data)]
  # data <- rename(data, colrenames)
  tryCatch(data <- rename(data, ...),
           error = function(e) stop("Error renaming data columns. Probably a column name does not exist.\nCheck whether an argument was misspelled and data name pairs are given in the form expected_name = true_name."))



  #### Check for opts given ####
  opts_missing <- !(c("nAttempts", "nRestarts", "maxfun", "maxit", "reltol", "factr") %in% names(opts))
  opts <- c(opts,
            setNames(as.list(c(5, 5, 8000, 2*1e+3,1e-6, 1e-10)[opts_missing]),
                     c("nAttempts", "nRestarts","maxfun", "maxit", "reltol", "factr")[opts_missing]))
  optim_method <- match.arg(optim_method, c("bobyqa", "L-BFGS-B", "Nelder-Mead"))

  if (!is.null(init_grid)) {
    if (nrow(init_grid)< opts$nAttempts) {
      opts$nAttempts <- nrow(init_grid)
      if (!opts_missing[1]) {
        message("opts$nAttempts is reduced to nrow(init_grid).")
      }
    }
  }
  if (!grid_search && is.null(init_grid) && !opts_missing[1]) {
    opts$nAttempts <- 1
    message("opts$nAttempts is reduced to 1, since grid_search was omitted.")
  }

  #### Check for column names given ####
  names_missing <- !(c("condition","response","stimulus","rating", "rt", "sbj", "correct") %in% names(data_names))
  data_names <- c(data_names,
                  setNames(as.list(c("condition","response","stimulus","rating", "rt", "sbj", "correct")[names_missing]),
                           c("condition","response","stimulus","rating", "rt", "sbj", "correct")[names_missing]))

  #### Process data input ####
  ### extract columns and check right formatting:
  cols <- names(data)
  condition <- data[[data_names$condition]]
  if(!is.factor(condition)) {
    condition <- factor(condition, levels=sort(unique(condition)))
  }
  if (data_names$stimulus %in% cols) {
    stimulus <- data[[data_names$stimulus]]
  }
  if (data_names$response %in% cols) {
    response <- data[[data_names$response]]
  }
  if (!exists("stimulus")) {
    if (!(data_names$correct %in% cols)) {stop("Column names in data must contain 2 of following 3: stimulus, response, correct or must be specified with the data_names argument.")}
    correct <- data[[data_names$correct]]
    stim_levels <- sort(unique(response))
    if (length(stim_levels) != 2) {stop(paste("Response must have exactly two unique values! Values are: ", paste(stim_levels, collapse = ", ")))}
    response <- if_else(response==stim_levels[1],-1,1)
    stimulus <- c(-1,1)[(1+as.numeric(response==1))]*(-1)^(1-correct)
  } else {
    if (!exists("response")) {
      if (!("correct" %in% cols)) {stop("Column names in data must contain 2 of following 3: stimulus, response, correct or must be specified with the data_names argument.")}
      correct <- data[[data_names$correct]]
      stim_levels <- sort(unique(stimulus))
      if (length(stim_levels) != 2) {stop(paste("Stimulus must have exactly two unique values! Values are: ", paste(stim_levels, collapse = ", ")))}
      stimulus <- if_else(stimulus==stim_levels[1],-1,1)
      response <- if_else(stimulus*(-1)^correct==1, -1, 1)
    } else {
      if (all(response %in% c(-1,1))) {
        if (!(all(stimulus %in% c(-1,1)))) {
          stim_levels <- sort(unique(stimulus))
          stimulus <- if_else(stimulus==stim_levels[1],-1,1)
        }
      } else {
        stim_levels <- sort(unique(response))
        response <- if_else(response==stim_levels[1],-1,1)
        if (!(all(stimulus %in% c(-1,1)))) {
          if (!(all(stimulus %in% stim_levels))) { stop("Values in stimulus must either be in c(-1,1) or same values as in response")}
          stimulus <- if_else(stimulus==stim_levels[1],-1,1)
        }
      }
    }
  }

  if (grepl("RM", model)) {         # If Race Models are used, re-code lower(-1) and upper(+1) boundary to first(1) and second(2) accumulator
    stimulus <- (stimulus/2 +1.5)   # Diffusion-based models (dynWEV; 2DSD): -1    1
    response <- (response/2 +1.5)   # Race-based models (IRM(t), PCRM(t)):    1    2
  }

  rt <- data[[data_names$rt]]
  rating <- data[[data_names$rating]]
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

  if (length(unique(rating))< nRatings) {
    # If some rating categories are not used, we fit less thresholds numerically and fill up the
    # rest by the obvious best-fitting thresholds (e.g. +/- Inf for the lowest/highest...) in the end.
    # Therefore, save the true rating-vector for later recovery of which thresholds are indeed fitted and
    # reduce the vector for the fitting process to a integer vector without gaps
    ### ToDo:   For sym_thetas==FALSE, use different nRatings for lower and upper responses in fitting
    used_cats <- sort(unique(rating))
    actual_nRatings <- nRatings
    rating <- as.integer(as.factor(rating))
    nRatings <- length(unique(rating))
  } else {
    used_cats <- NULL
    actual_nRatings <- NULL
  }
  if ( nRatings < 2) {
    stop("There has to be at least two rating levels")
  }
  nConds <- length(levels(condition))
  if (nConds ==0) {
    condition = 1
    nConds = 1
  }
  df = data.frame(rating, stimulus, response, condition,rt)
  df$n = 1
  df = aggregate(n~rating+stimulus+response+condition+rt, df, sum)


  #### Initialize logger, if logging is wished ####
  if (logging==TRUE) {
    ## get participant ID for logging
    if (data_names$sbj %in% cols){
      sbjcol <- data_names$sbj
    } else if ("participant" %in% cols) {
      sbjcol <- "participant"
    }
    if (!exists("sbjcol")) {
      participant <- 999
    } else {
      if (length(unique(c(t(data[,sbjcol]))))==1) {
        participant <- as.numeric(data[1,sbjcol])
      } else {
        participant <- 999
      }
    }
    ## set logger and logging file
    dir.create("autosave", showWarnings = FALSE)
    dir.create(paste("autosave/fit", model, sep=""), showWarnings = FALSE)
    filename = paste("autosave/fit", model,"/part_", participant,".RDATA", sep = "")
    logger <- layout_glue_generator(format = paste('{level} [{time}] on process {pid} {fn} for participant ',participant,' and model ', model,': {msg}', sep=""))
    log_layout(logger)
    log_appender(appender_file(file=paste("autosave/fit", model,"/logging_", model, ".txt", sep="")), index=2)
    log_threshold(DEBUG, index=2)
    log_threshold(DEBUG)
    log_layout(logger, index=2)
  }


  ### Adapt arguments for individual settings
  if (!is.list(fixed)) fixed <- as.list(fixed)
  if (!("sym_thetas" %in% names(fixed))) fixed["sym_thetas"] <- FALSE
  sym_thetas <- as.logical(fixed['sym_thetas'])
  fixed <- fixed[names(fixed)!="sym_thetas"]

  ### Now, call the specific fitting functions:
  if (model == "2DSD") res <- fitting2DSD(df, nConds, nRatings, fixed, sym_thetas,
                                          grid_search, init_grid, optim_method, opts,
                                          logging, filename,
                                          useparallel, n.cores,
                                          restr_tau, precision,
                                          used_cats, actual_nRatings)
  if (model == "dynWEV") res <- fittingdynWEV(df, nConds, nRatings, fixed, sym_thetas,
                                              grid_search, init_grid, optim_method, opts,
                                              logging, filename,
                                              useparallel, n.cores,
                                              restr_tau, precision,
                                              used_cats, actual_nRatings)
  if (grepl("IRM", model)) res <- fittingIRM(df, nConds, nRatings, fixed,
                                             sym_thetas, grepl("t", model),
                                        grid_search, init_grid, optim_method, opts,
                                        logging, filename,
                                        useparallel, n.cores,
                                        used_cats, actual_nRatings)
  if (grepl("PCRM", model)) res <- fittingPCRM(df, nConds, nRatings, fixed ,
                                               sym_thetas, grepl("t", model),
                                          grid_search, init_grid, optim_method, opts,
                                          logging, filename,
                                          useparallel, n.cores,
                                          used_cats, actual_nRatings)

  return(res)
}
