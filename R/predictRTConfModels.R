#' Prediction of confidence and RT distributions for several sequential
#' sampling confidence models and parameter constellations in parallel
#'
#' This function is a wrapper around the functions \code{\link{predictRTConf}} (see
#' there for more information). It calls the respective function for predicting the
#' response distribution (discrete decision and rating outcomes) and the rt density
#' (density for decision, rating and response time) for every model and
#' participant combination in \code{paramDf}.
#' Also, see \code{\link{dWEV}}, \code{\link{d2DSD}}, and \code{\link{dRM}} for more
#' information about the parameters.
#'
#' @param paramDf a dataframe with one row per combination of model and
#' participant/parameter set. Columns must include a participant or subject column,
#' a model column and the names of the model parameters. For different stimulus
#' quality/mean drift rates, names should be v1, v2, v3,.... Different s parameters
#' are possible with s1, s2, s3... with equally many steps as for drift rates (same
#' for sv parameter in dynWEV and 2DSD).
#' Additionally, the confidence thresholds should be given by names with
#' thetaUpper1, thetaUpper2,..., thetaLower1,... or,
#' for symmetric thresholds only by theta1, theta2,....
#' @param maxrt numeric. The maximum RT for the
#' integration/density computation. Default: 15 (for \code{predictConfModels} (integration)) and
#' 9 (for \code{predictRTModels}).
#' @param subdivisions \code{integer} (default: 100).
#' For \code{predictConfModels} it is used as argument for the inner integral routine.
#' For \code{predictRTModels} it is the number of points for which the density is computed.
#' @param minrt numeric or NULL(default). The minimum rt for the density computation.
#' If NULL, the minimal possible response time possible with given parameters will be used (min(t0)).
#' @param  simult_conf logical, only relevant for dynWEV and 2DSD. Whether in the experiment
#' confidence was reported simultaneously with the decision, as then decision and confidence
#' judgment are assumed to have happened subsequent before response and computations are
#' different, when there is an observable interjudgment time (then `simult_conf` should be FALSE).
#' @param scaled logical. Whether the computed density should be scaled to integrate to one
#' (additional column densscaled). Otherwise the output is a defective density (i.e. its
#' integral is equal to the probability of a response and not 1). If TRUE, the argument
#' `DistConf` should be given, if available. Default: FALSE.
#' @param DistConf NULL or data.frame. A data.frame with participant
#' and model columns and columns, giving the distribution of response and rating choices for
#' different conditions and stimulus categories in the form of the output of
#' \code{predictConfModels}. It is only necessary if `scaled=TRUE`, because these
#' probabilities are used for scaling. If `scaled=TRUE` and `DistConf=NULL`, it will be computed
#' with the function \code{predictConfModels}, which takes some time and the function will
#' throw a message. Default: NULL
#' @param stop.on.error logical. Argument directly passed on to integrate. Default is FALSE,
#' since the densities invoked may lead to slow convergence of the integrals (which are still
#' quite accurate) which causes R to throw an error.
#' @param parallel logical. If TRUE, prediction is parallelized over participants and models
#' (i.e. over the calls for the respective \code{\link{predictRTConf}} functions).
#' @param n.cores integer. If \code{parallel} is TRUE, the number of cores used for
#' parallelization is required. If NULL (default) the number of available cores -1 is used.
#' @param .progress logical. If TRUE (default) a progress bar is drawn to the console. (Works
#' for some OS only when `parallel=FALSE`.)
#'
#' @return \code{predictConfModels} gives a data frame/tibble with columns: participant, model,
#' condition, stimulus, response, rating, correct, p, info, err. p is the predicted probability
#' of a response and rating, given the stimulus category and condition. Message and error refer
#' to the respective outputs of the integration routine used for computation.
#' \code{predictRTModels} returns a data frame/tibble with columns: participant, model,
#' condition, stimulus, response, rating, correct, rt and dens (and densscaled, if `scaled=TRUE`).
#'
#'
#' @details These functions merely split the input data frame by model participants combinations,
#' call the equivalent \code{\link{predictRTConf}} functions for the individual parameter sets
#' and bind the outputs together. They are included for convenience and the easy parallelization,
#' which facilitates speeding up computations considerably. For the argument
#' \code{paramDf}, the output of the fitting function \code{\link{fitRTConfModels}} with the
#' respective models and participants may be used.
#'
#' The function \code{\link{predictConf}} (called by \code{predictConfModels})
#' consists merely of an integration of the reaction time density or the given model,
#' \code{{d*model*}}, over the reaction time in a reasonable interval (0 to maxrt).
#' The function \code{\link{predictRT}} (called by \code{predictRTModels}) wraps these
#' density functions to a parameter set input and a data.frame output. '
#' Note, that the encoding for stimulus identity is different between diffusion based models
#' (2DSD, dynWEV) and race models (IRM(t), PCRM(t)). Therefore, in the columns stimulus and
#' response there will be a mix of encodings: -1/1 for diffusion based models and 1/2 for
#' race models. This, usually is not important, since for further aggregation models will
#' not be mixed.
#'
#' @note Different parameters for different conditions are only allowed for drift rate
#' \code{v}, drift rate variability \code{sv} (only dynWEV and 2DSD), and process variability
#' \code{s}. All other parameters are used for all conditions.
#'
#'
#' @author Sebastian Hellmann.
#'
#' @name predictRTConfModels
#' @importFrom dplyr rename
#' @import parallel
#' @importFrom parallel detectCores
# @importFrom pracma integral
#' @aliases predictConfModels
#' @importFrom Rcpp evalCpp
#'

#' @rdname predictRTConfModels
#' @export
predictConfModels <- function(paramDf,
                              maxrt=15, subdivisions = 100L, simult_conf = FALSE,
                              stop.on.error=FALSE,
                              .progress=TRUE,
                              parallel = TRUE, n.cores=NULL){ #  ?ToDO: vary_sv=FALSE, RRT=NULL, vary_tau=FALSE
  models <- unique(paramDf$model)
  if (!is.numeric(maxrt)) stop("maxrt must be numeric")
  if (!all(models %in% c("IRM", "PCRM", "IRMt", "PCRMt", "dynWEV", "2DSD"))) stop("model must be 'dynWEV', '2DSD', 'IRM', 'PCRM', 'IRMt', or 'PCRMt'")
  sbjcol <- c("subject", "participant", "sbj")[which(c("subject", "participant", "sbj") %in% names(paramDf))]
  if (length(sbjcol)==0) {
    paramDf$sbj <- 1
    sbjcol <- "sbj"
  } else {
    if (sbjcol != "sbj") {
      paramDf$sbj <- paramDf[[sbjcol]]
      paramDf[[sbjcol]] <- NULL
    }
  }
  subjects <- unique(paramDf$sbj)


  ### Determine number of jobs, i.e. model-participant-combinations
  nJobs <- nrow(paramDf)

  call_predfct <- function(X) {
    cur_model <- models[X[1]]
    cur_sbj <- X[2]
    sbj <- NULL # to omit a note because of an unbound variable
    model <- NULL # to omit a note because of an unbound variable
    params <- subset(paramDf, sbj==cur_sbj & model==cur_model)
    res <- predictConf(params, cur_model,
                       maxrt=maxrt, subdivisions = subdivisions,
                       simult_conf = simult_conf,
                       stop.on.error=stop.on.error,
                       .progress=.progress)
    res$model <- cur_model
    res[[sbjcol]] <- cur_sbj
    return(res)
  }

  jobs <- expand.grid(model=1:length(models), sbj=subjects)

  if (parallel) {
    listjobs <- list()
    for (i in 1:nrow(jobs)) {
      listjobs[[i]] <- c(model = jobs[["model"]][i], sbj = jobs[["sbj"]][i])
    }
    if (is.null(n.cores)) {
      n.cores <- detectCores()-1
    }
    n.cores <- min(n.cores, nJobs)
    cl <- makeCluster(type="SOCK", n.cores)
    clusterExport(cl, c("paramDf", "sbjcol", "models",
                              "subjects", "maxrt",
                              "subdivisions", "simult_conf", "stop.on.error",
                              ".progress", "call_predfct"), envir = environment())
    on.exit(try(stopCluster(cl), silent = TRUE))
    res <- clusterApplyLB(cl, listjobs, fun=call_predfct)
    stopCluster(cl)
  } else {
    res <- apply(jobs, 1, call_predfct)
  }
  res <- do.call(rbind, res)
  return(res)
}






#' @rdname predictRTConfModels
#' @export
predictRTModels <- function(paramDf,
                            maxrt=9, subdivisions = 100L,  minrt=NULL,
                            simult_conf = FALSE,
                            scaled = FALSE, DistConf=NULL,
                            .progress = TRUE,
                            parallel = TRUE, n.cores=NULL){ #  ?ToDO: vary_sv=FALSE, RRT=NULL, vary_tau=FALSE
  if (!is.numeric(maxrt)) stop("maxrt must be numeric")
  models <- unique(paramDf$model)
  if (!all(models %in% c("IRM", "PCRM", "IRMt", "PCRMt", "dynWEV", "2DSD"))) stop("model must be 'dynWEV', '2DSD', 'IRM', 'PCRM', 'IRMt', or 'PCRMt'")
  sbjcol <- c("subject", "participant", "sbj")[which(c("subject", "participant", "sbj") %in% names(paramDf))]
  if (length(sbjcol)==0) {
    paramDf$sbj <- 1
    sbjcol <- "sbj"
  } else {
    if (sbjcol != "sbj") {
      paramDf$sbj <- paramDf[[sbjcol]]
      paramDf[[sbjcol]] <- NULL
    }
  }
  subjects <- unique(paramDf$sbj)

  if (scaled && !is.null(DistConf)) {
    if (sbjcol != "sbj") {
      DistConf$sbj <- DistConf[[sbjcol]]
      DistConf[[sbjcol]] <- NULL
    }
    if (!all(subjects %in% unique(DistConf$sbj))) stop("There is a subject in paramDf, that is not present in DistConf!")
    if (!all(models %in% unique(DistConf[["model"]]))) stop("There is a model in paramDf, that is not present in DistConf!")
  }

  if (scaled && is.null(DistConf)) {
    message(paste0("scaled is TRUE and DistConf is NULL.\n",
                   "Confidence distribution is calculated before ",
            "computing the RT densities,\nthis takes considerable additional time..."))
    DistConf <- predictConfModels(paramDf, maxrt = 15,
                                  subdivisions = 100L,
                                  stop.on.error=FALSE,
                                  simult_conf = simult_conf,
                                  .progress=.progress,
                                  parallel = parallel, n.cores=n.cores)
    #if (is.null(return_DistConf)) return_DistConf <- TRUE
    message("...finished computation of confidence distribution.")
  }



  ### Determine number of jobs, i.e. model-participant-combinations
  nJobs <- nrow(paramDf)
  ### set up parallelization settings
  if (parallel) {
    if (is.null(n.cores)) {
      n.cores <- detectCores()-1
    }
    n.cores <- min(n.cores, nJobs)
  }
  if (is.null(minrt)) {
    minrt <- min(paramDf$t0)
    if (simult_conf && any(c("2DSD", "dynWEV") %in% models)) {
      pars_diffmodels <- filter(paramDf, .data$model %in% c("2DSD", "dynWEV"))
      minrt <- min(minrt, pars_diffmodels$t0+pars_diffmodels$tau)
    }
  }
  call_predfct <- function(X) {
    cur_model <- models[X[1]]
    cur_sbj <- X[2]
    sbj <- NULL # to omit a note because of an unbound variable
    model <- NULL # to omit a note because of an unbound variable
    params <- subset(paramDf, sbj==cur_sbj & model==cur_model)
    if (scaled) {
      cur_DistConf <- subset(DistConf, sbj==cur_sbj & model==cur_model)
    } else {
      cur_DistConf <- NULL
    }
    res <- predictRT(params, cur_model,
                     maxrt=maxrt, subdivisions = subdivisions,
                     minrt=minrt, simult_conf = simult_conf,
                     scaled=scaled, DistConf=cur_DistConf,
                     .progress=.progress)
    res$model <- cur_model
    res[[sbjcol]] <- cur_sbj
    return(res)
  }
  jobs <- expand.grid(model=1:length(models), sbj=subjects)

  if (parallel) {
    listjobs <- list()
    for (i in 1:nrow(jobs)) {
      listjobs[[i]] <- c(model = jobs[["model"]][i], sbj = jobs[["sbj"]][i])
    }
    cl <- makeCluster(type="SOCK", n.cores)
    clusterExport(cl, c("paramDf", "sbjcol", "models",
                              "subjects", "maxrt", "minrt", "simult_conf",
                              "subdivisions", "scaled", "DistConf",
                              ".progress", "call_predfct"), envir = environment())
    on.exit(try(stopCluster(cl), silent = TRUE))
    res <- clusterApplyLB(cl, listjobs, fun=call_predfct)
    stopCluster(cl)
  } else {
    res <- apply(jobs, 1, call_predfct)
  }
  res <- do.call(rbind, res)
  return(res)
}