#' Function for fitting sequential sampling confidence models with custom manipulations
#'
#' Fits the parameters of different models of response time and confidence, including
#' the 2DSD model (Pleskac & Busemeyer, 2010), dynWEV, DDMConf, and various
#' flavors of race models (Hellmann et al., 2023), where manipulations are
#' specified using formula for each manipulated parameter.
#' Which model to fit is
#' specified by the argument \code{model}.
#' Only a ML method is implemented.
#' See \code{\link{dWEV}}, \code{\link{d2DSD}}, and \code{\link{dRM}} for more
#' information about the parameters and Details for not-fitted parameters.
#'
#' @param data a `data.frame` where each row is one trial, containing at
#' least the following variables (others may be used in formulae):
#' * \code{rating} (discrete confidence judgments, should be given as integer vector; otherwise will be transformed to integer),
#' * \code{rt} (giving the reaction times for the decision task),
#' * either 2 of the following (see details for more information about the accepted formats):
#'   * \code{stimulus} (encoding the stimulus category in a binary choice task),
#'   * \code{response} (encoding the decision response),
#'   * \code{correct} (encoding whether the decision was correct; values in 0, 1)
#' * \code{sbj} or \code{participant} (optional; giving the subject ID; only relevant if `logging = TRUE`;
#'                                                       if unique the ID is used in saved files with interim results
#'                                                       and logging messages;
#'                                                       if non-unique or missing and `logging =TRUE`, 999 will be used then)
#' @param model character scalar. One of "dynWEV", "2DSD", "IRM", "PCRM", "IRMt", "PCRMt", or "DDMConf" for the model to be fit.
#' @param fixed list. List with parameter-value pairs for parameters that should not be fitted. See Details.
#' @param init_grid data.frame or `NULL`. Grid for the initial parameter search. Each row is one parameter constellation.
#' See details for more information. If \code{NULL} a default grid will be used.
#' @param grid_search logical. If `FALSE`, the grid search before the optimization
#' algorithm is omitted. The fitting is then started with a mean parameter set
#' from the default grid (if `init_grid=NULL`) or directly with the rows from
#' `init_grid`, if not `NULL`. (Default: `TRUE`)
#' @param nRatings integer. Number of rating categories. If `NULL`, the maximum of
#' `rating` and `length(unique(rating))` is used. This argument is especially
#' important for data sets where not the whole range of rating categories is realized.
#' If given, ratings has to be given as factor or integer.
#' @param logging logical. If `TRUE`, a folder 'autosave/fit**model**' is created and
#' messages about the process are printed in a logging file and to console (depending
#' on OS). Additionally intermediate results are saved in a `.RData` file with the
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
#' @param restr_tau numerical or `Inf` or `"simult_conf"`. For 2DSD and dynWEV only.
#' Upper bound for tau. Fits will be in the interval (0,`restr_tau`). If FALSE tau will be unbound.
#' For `"simult_conf"`, see the documentation of \code{\link{d2DSD}} and \code{\link{dWEV}}
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
#'  Correct must always be given as a 0-1-vector. If the stimulus column is given
#'  together with a response column and they both do not match the above format,
#'  they need to have the same values/levels (if `factor`).
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
#'  \code{z=.5} in a diffusion process.
#'
#'  \strong{Parameters not fitted}. The models get developed continuously and not
#'  all changes are adopted in the fitting function instantly. Following parameters
#'  are currently not included in the fitting routine:
#'  - in race models: \code{sza}, \code{szb}, \code{smu1}, and \code{smu2}
#'
#' \strong{`init_grid`}. Each row should be one parameter set to check. The column names
#' should include the parameters of the desired model, which are the following for 2DSD:
#' `a`, `vmin` and `vmax` (will be equidistantly spanned across conditions), `sv`, `z` (as the
#' relative starting point between 0 and `a`), `sz` (also in relative terms), `t0`, `st0`, `theta0`
#' (minimal threshold), `thetamax` (maximal threshold; the others will be equidistantly
#' spanned symmetrically for both decisions), and `tau`. For dynWEV,
#' additionally `w` , `svis`, and `sigvis` are required. For the race models the parameters
#' are: `vmin`, `vmax` (will be equidistantly
#' spanned across conditions), `a` and `b` (decision thresholds), `t0`, `st0`, `theta0`
#' (minimal threshold), `thetamax` (maximal threshold;
#' the others will be equidistantly spanned symmetrically for both decisions), and for
#' time-dependent confidence race models
#' additionally `wrt` and `wint` (as weights compared to `wx=1`).
#'
#'  \strong{opts}. A list with numerical values. Possible options are listed below
#'  (together with the optimization method they are used for).
#'  * \code{nAttempts} (all) number of best performing initial parameter sets used for
#'   optimization; default 5, if `grid_search` is `TRUE`.
#'  If `grid_search` is `FALSE` and `init_grid` is `NULL`, then `nAttempts` will be set to 1 (and
#'  any input will be ignored).
#'  If `grid_search` is `FALSE` and `init_grid` is not `NULL`, the rows of `init_grid` will be used
#'  from top to bottom
#'  (since no initial grid search is done) with not more than `nAttempts` rows used.
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
#' <https://nashjc.wordpress.com/2016/11/10/why-optim-is-out-of-date/>
#'
#' <https://www.damtp.cam.ac.uk/user/na/NA_papers/NA2009_06.pdf>
#'
#'
#'
#' @author Sebastian Hellmann.
#'
#' @name fitRTConf_formula
#' @importFrom stats setNames aggregate optim qnorm pnorm optimize quantile
#' @importFrom minqa bobyqa
#' @importFrom dplyr if_else rename
#' @importFrom rlang .data
#' @importFrom magrittr %>%
#' @import parallel
# @importFrom pracma integral
#' @aliases fitSeqSampConf fitConfModel fitConf fitConfRT
#'
#' @examples
#' # We use one of the implemented models, "dynWEV"
#' # 1. Generate data
#' # data with positive drift (stimulus = "upper")
#' data <- rWEV(20, a=2,v=0.5,t0=0.2,z=0.5, sz=0.1,sv=0.1, st0=0,  tau=4, s=1, w=0.3)
#' data$stimulus <- "upper"
#' # data with negtive drift (stimulus = "lower") but same intensity
#' data2 <- rWEV(100, a=2,v=-0.5,t0=0.2,z=0.5,sz=0.1,sv=0.1, st0=0,  tau=4, s=1, w=0.3)
#' data2$stimulus <- "lower"
#' data <- rbind(data, data2)
#' # Transfer response column and add dummy condition column
#' data$response <- ifelse(data$response==1, "upper", "lower")
#' data$condition <- 1
#' # Take some confidence thresholds for discrete ratings
#' threshs <- c(-Inf, 1, 2, Inf)
#' data$rating <- as.numeric(cut(data$conf, breaks = threshs, include.lowest = TRUE))
#' head(data)
#'
#' # 2. Use fitting function
#' # Fitting the model with these opts results in a pretty bad fit
#' # (especially because of omitting the grid_search)
#' \donttest{
#'    fitRTConf(data, "dynWEV", fixed=list(sym_thetas=TRUE, z=0.5, st0=0),
#'             grid_search = FALSE, logging=FALSE,
#'             opts = list(nAttempts=1, nRestarts=2, maxfun=2000))
#'  }
#'



#' @rdname fitRTConf_formula
#' @export
fitRTConf_formula <- function(data, model = "dynWEV",
                              manipulations = list(),
                              fixed = list(sym_thetas = FALSE, s=1),
                              nRatings = NULL, restr_tau =Inf,
                              precision=1e-5,grid_search=TRUE,
                              logging=FALSE, opts=list(), optim_method = "Nelder-Mead",
                              useparallel = FALSE, n.cores=NULL, ...){ #  ?ToDO: vary_sv=FALSE, RRT=NULL, vary_tau=FALSE
  # Check if package 'logger' is installed, if logging is wished
  if (logging && !requireNamespace("logger", quietly = TRUE)) {
    warning("Package 'logger' is not installed but needed to log fitting progress.
            Process continues withouth logging.
            Interrupt and install 'logger' if logging is needed.", immediate.=TRUE)
    logging <- FALSE
  }

  #### Check for opts given ####
  opts_missing <- !(c("nAttempts", "nRestarts", "maxfun", "maxit", "reltol", "factr") %in% names(opts))
  opts <- c(opts,
            setNames(as.list(c(5, 5, 8000, 2*1e+3,1e-6, 1e-10)[opts_missing]),
                     c("nAttempts", "nRestarts","maxfun", "maxit", "reltol", "factr")[opts_missing]))
  optim_method <- match.arg(optim_method, c("bobyqa", "L-BFGS-B", "Nelder-Mead"))

  #### Process data input ####
  ### extract columns and check right formatting:
  cols <- names(data)
  if ("stimulus" %in% cols) {
    stimulus <- data[["stimulus"]]
  }
  if ("response" %in% cols) {
    response <- data[["response"]]
  }
  if (!exists("response")) {
    if (!("correct" %in% cols)) {stop("Column names in data must contain 'response' or 'correct'.")}
    if (!exists("stimulus")) {
      if (!("z" %in% fixed) || !(fixed$sym_thetas)) warning("Only the accuracy of responses provided. We recommend to set z=0.5 and sym_thetas=TRUE in the 'fixed' argument in such situations because biases could not be identified!")
      response <- data[["correct"]]
    } else {
      stim_levels <- sort(unique(stimulus))
      if (length(stim_levels) != 2) {stop(paste("Stimulus must have exactly two unique values! Values are: ", paste(stim_levels, collapse = ", ")))}
      stimulus <- if_else(stimulus==stim_levels[1],-1,1)
      response <- if_else(stimulus*(-1)^correct==1, -1, 1)
    }
  }
  if (all(response %in% c(0,1))) {
    response <- as.logical(response)
  } else {
    response <- response == -1
  }

  rt <- data[["rt"]]
  rating <- data[["rating"]]
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
    if (max(length(unique(rating)), max(rating)+1-min(rating))==nRatings && any(rating==0)) {
      rating <- rating + 1L
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
  DVs = data.frame(rating,response,rt)

  ### Use manipulations and fixed parameters to get a model.matrix for fitting

  # define all necessary parameters for a model
  model_name <- model
  models <- c("2DSD", "2DSDT", "dynWEV", "dynaViTE")
  matches <- sapply(models, function(m) grepl(m, model))
  model <- models[max(which(matches))]
  parnames <- c("v", "z", "a", "d", "sz", "t0", "st0", "sv", "tau", "w", "svis", "sigvis", "lambda", "s")
  model_params_fixed <- c("lambda"=0, "w"=1, "sigvis"=1, "svis"=1)
  if (model=="2DSD") {
    to_set <- setdiff(c("w", "sigvis", "svis", "lambda"), names(fixed))
  } else if (model=="2DSDT") {
    to_set <- setdiff(c("w", "sigvis", "svis"), names(fixed))
  } else if (model=="dynWEV") {
    to_set <- setdiff(c("lambda"), names(fixed))
  } else if (model=="dynaViTE") {
    to_set <- NULL
  } else { stop(paste0("Model: ", model, " not implemented!")) }
  fixed <- c(fixed, model_params_fixed[to_set])

  # Get manipulated parameter names
  man_pars <- sapply(manipulations, FUN = function(x) x[[2]])
  const_pars <- setdiff(parnames, c(man_pars, names(fixed)))
  names(manipulations) <- unlist(man_pars)

  # # Create one vector with a formula for each model parameter,
  # # those that are not manipulated get only an intercept
  # # if a parameter is fixed, it gets no formula
  # manipulations2 <- list()
  # for (i in 1:length(parnames)) {
  #   # Only, if parameter is not excluded from fitting (i.e. fixed)
  #   if (!parnames[i] %in%names(fixed)) {
  #     # If parameter is manipulated, use user-specified formula
  #     if (parnames[i] %in% man_pars) {
  #       manipulations2[[parnames[i]]] <- manipulations[[match(parnames[i], man_pars)]]
  #     } else {
  #       # else use only an intercept
  #       manipulations2[[parnames[i]]] <- as.formula(paste0(parnames[i], "~1"))
  #     }
  #     # create a constant column in the data for each parameter (necessary for model.matrix to function)
  #     data[[parnames[i]]] <- 1
  #   } else {
  #     if (parnames[i] %in% man_pars) warning(paste0("Parameter `", parnames[i], "` is ambiguously specified as both fixed and manipulated!\nIt will be fixed for fitting!"))
  #   }
  # }
  # # drop Input and further use filled manipulation vector
  # manipulations <- manipulations2

  # create one big model matrix for all manipulation
  # containing all necessary columns for any factor level etc.
  model_matrix <- matrix(nrow = nrow(data), ncol=0)
  # furhter, create a list where each parameter gets an element,
  # containing a character vector with the column names of the
  # model.matrix used for this parameter
  fit_pars_columns <- list()
  fit_pars_beta_names <- NULL
  fit_par_names <- NULL
  # create a constant column in the data for each parameter (necessary for model.matrix to function)
  data[,names(manipulations)] <- 1

  factor_cols <- names(data[,sapply(data, is.factor), drop=FALSE])
  #contr_list <- eval(parse(text=paste0("list(", paste(paste(factor_cols, "='contr.treatment'"), c(rep(",", length(factor_cols)-1), " "), collapse = ""), ")")))
  contr_list <- NULL

  for (i in 1:length(manipulations)) {
    cur_parameter <- as.character(manipulations[[i]][[2]])
    #model.frame(manipulations[[i]], data=data)
    temp <- model.matrix(manipulations[[i]],data,
                         contrasts.arg = contr_list )#, contrasts.arg =
    fit_pars_columns[[names(manipulations)[i]]] <- colnames(temp)
    fit_pars_beta_names <- c(fit_pars_beta_names, paste(cur_parameter, colnames(temp), sep="_"))
    fit_par_names <- c(fit_par_names, cur_parameter)
    model_matrix <- cbind(model_matrix, temp[, setdiff(colnames(temp), colnames(model_matrix)), drop=FALSE])
  }
  # model_matrix
  # names(fit_pars_columns) <- fit_par_names
  # fit_pars_columns
  # #lapply(seq_along(fit_pars_columns), function(y, n, i) { paste(n[[i]], y[[i]], sep="_") }, y=fit_pars_columns, n=names(fit_pars_columns))
  # fit_par_names
  # #fit_pars_columns_names <- sub("\\(Intercept\\)", "1",fit_pars_columns_names)
  # fit_pars_beta_names
  # const_pars



  #### Initialize logger, if logging is wished ####
  if (logging==TRUE) {
    ## get participant ID for logging
    if ("sbj" %in% cols){
      sbjcol <- "sbj"
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
    dir.create(paste("autosave/fit", model_name, sep=""), showWarnings = FALSE)
    filename = paste("autosave/fit", model_name,"/part_", participant,".RDATA", sep = "")
    logger <- logger::layout_glue_generator(format = paste('{level} [{time}] on process {pid} {fn} for participant ',participant,' and model ', model_name,': {msg}', sep=""))
    logger::log_layout(logger)
    logger::log_appender(logger::appender_file(file=paste("autosave/fit", model_name,"/logging_", model_name, ".txt", sep="")), index=2)
    logger::log_threshold(logger::DEBUG, index=2)
    logger::log_threshold(logger::DEBUG)
    logger::log_layout(logger, index=2)
  }


  ### Adapt arguments for individual settings
  if (!is.list(fixed)) fixed <- as.list(fixed)
  if (!("sym_thetas" %in% names(fixed))) fixed["sym_thetas"] <- FALSE
  sym_thetas <- as.logical(fixed['sym_thetas'])
  fixed <- fixed[names(fixed)!="sym_thetas"]

  ### Now, call the specific fitting functions:
  if (grepl("2DSD", model)) {
    if (model=="2DSD") fixed$lambda <- 0
    res <- fitting2DSD(df, nConds, nRatings, fixed, sym_thetas,
                       grid_search, init_grid, optim_method, opts,
                       logging, filename,
                       useparallel, n.cores,
                       restr_tau, precision,
                       used_cats, actual_nRatings)
  }
  if (grepl("dynWEV|dynaViTE",model)) {
    if (model=="dynWEV") fixed$lambda <- 0
    res <- fittingdynWEV_formula(DVs, model_matrix,
                                 fixed,
                                 fit_pars_beta_names, const_pars,
                                 fit_pars_columns,
                                 nRatings, sym_thetas,
                                 optim_method, opts,
                                 logging, filename,
                                 useparallel, n.cores,
                                 restr_tau, precision, grid_search,
                                 used_cats, actual_nRatings)
  }
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
  if (model == "DDMConf") res <- fittingDDMConf(df, nConds, nRatings, fixed, sym_thetas,
                                                grid_search, init_grid, opts,
                                                logging, filename,
                                                useparallel, n.cores,
                                                precision,
                                                used_cats, actual_nRatings, precision)
  if (!exists("res")) stop("model not known. model must contain one of: 'dynaViTE', 'dynWEV', '2DSD', '2DSDT', 'IRM', 'PCRM', or 'DDMConf'")


  return(res)
}
