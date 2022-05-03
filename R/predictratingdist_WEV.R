#' Prediction of Rating and Reaction Time Distribution in dWEV and 2DSD confidence models
#'
#' \code{predictWEV_Conf} predicts the categorical response distribution of
#' decision and confidence ratings, \code{predictWEV_RT} computes the predicted
#' RT-Distribution (density) in the 2DSD Model (Pleskac & Busemeyer, 2010) and the
#' dWEV model (Hellmann & Rausch), given specific parameter constellations.
#' See \code{\link{dWEV}} and \code{\link{d2DSD}} for more information about parameters.
#'
#' @param paramDf a list or dataframe with one row. Column names should match the
#' names of dWEV and 2DSD model specific parameter names. For different stimulus
#' quality/mean drift rates, names should be v1, v2, v3,.... Different sv parameters
#' are possible with sv1, sv2, sv3... with equally many steps as for drift rates.
#' Additionally, the confidence thresholds should be given by names with
#' thetaUpper1, thetaUpper2,..., thetaLower1,... or,
#' for symmetric thresholds only by theta1, theta2,....
#' @param model character scalar. One of "dynWEV" or "2DSD".
#'
#' @param precision \code{numerical} scalar value. Precision of calculation.
#' Corresponds roughly to the number of decimals of the predicted CDFs that are
#' calculated accurately. Default is 3.
#' @param maxrt numeric. The maximum RT for the
#' integration/density computation. Default: 15 (for _Conf (integration)), 9 (for _RT).
#' @param subdivisions \code{integer} (default: 100).
#' For \code{predictWEV_Conf} it is used as argument for the inner integral routine.
#' For \code{predictWEV_RT} it is the number of points for which the density is computed.
#' @param minrt numeric or NULL(default). The minimum rt for the density computation.
#' @param  simult_conf logical. Whether in the experiment confidence was reported simultaneously
#' with the decision, as then decision and confidence judgment are assumed to have happened
#' subsequent before response and computations are different, when there is an observable
#' interjudgment time (then simult_conf should be FALSE).
#' @param scaled logical. For \code{predictWEV_RT}. Whether the computed density
#' should be scaled to integrate to one (additional column densscaled). Otherwise the output
#' is a defective density (i.e. its integral is equal to the probability of a response and
#' not 1). If TRUE, the argument DistConf should be given, if available. Default: FALSE.
#' @param DistConf NULL or data.frame. For \code{predictWEV_RT}. A data.frame or matrix
#' with column names, giving the distribution of response and rating choices for
#' different conditions and stimulus categories in the form of the output of
#' \code{predictWEV_Conf}. It is only necessary, if scaled=TRUE, because these
#' probabilities are used for scaling. If scaled=TRUE and DistConf=NULL, it will be computed
#' with the function \code{predictWEV_Conf}, which takes some time and the function will
#' throw a message. Default: NULL
#' @param stop.on.error logical. Argument directly passed on to integrate. Default is FALSE, since
#' the densities invoked may converge slowly (but are still quite accurate) which causes R to throw an
#' error.
#' @param .progress logical. if TRUE (default) a progress bar is drawn to the console.
#'
#' @return \code{predictWEV_Conf} gives a data frame/tibble with columns: condition, stimulus,
#' response, rating, correct, p, info, err. p is the predicted probability of a response
#' and rating, given the stimulus category and condition. Message and error refer to the
#' respective outputs of the integration routine used for prediction.
#' \code{predictWEV_RT} returns a data frame/tibble with columns: condition, stimulus,
#' response, rating, correct, rt and dens (and densscaled, if scaled=TRUE).
#'
#'
#' @details The function \code{predictWEV_Conf} consists merely of an integration of
#' the reaction time density, \code{\link{dWEV}} and \code{\link{d2DSD}}, over the reaction time in a reasonable
#' interval (\code{t0} to maxrt). The function \code{predictWEV_RT} wraps these density
#' functions to a parameter set input and a data.frame output.
#' For the argument \code{paramDf}, the output of the fitting function \code{\link{fitRTConf}}
#' with the respective model may be used.
#'
#' @note Different parameters for different conditions are only allowed for drift rate,
#' \code{v}, and drift rate variability, \code{sv}. All other parameters are used for all
#' conditions.
#'
#' @references Pleskac, T. J., & Busemeyer, J. R. (2010). Two-Stage Dynamic Signal Detection:
#' A Theory of Choice, Decision Time, and Confidence, \emph{Psychological Review}, 117(3),
#' 864-901. doi:10.1037/a0019737
#'
#' Rausch, M., Hellmann, S., & Zehetleitner, M. (2018). Confidence in masked orientation
#' judgments is informed by both evidence and visibility. \emph{Attention, Perception, &
#' Psychophysics}, 80(1), 134–154.  doi: 10.3758/s13414-017-1431-5
#'
#'
#' @author Sebastian Hellmann.
#'
#' @name predictWEV
#' @importFrom stats integrate
#' @import dplyr
#' @importFrom progress progress_bar
#' @importFrom magrittr %>%
#' @importFrom rlang .data
# @importFrom pracma integral
#' @aliases predictWEV predictdynWEV predict2DSD
#' @importFrom Rcpp evalCpp
#'

#' @rdname predictWEV
#' @export
predictWEV_Conf <- function(paramDf, model="dynWEV",
                                precision=1e-5,
                                maxrt=15, subdivisions = 100L, simult_conf = FALSE,
                            stop.on.error = FALSE,
                                .progress=TRUE){
  if (model =="WEVmu") model <- "dynWEV"
  nConds <- length(grep(pattern = "^v[0-9]", names(paramDf), value = T))
  symmetric_confidence_thresholds <- length(grep(pattern = "thetaUpper", names(paramDf), value = T))<1
  if (symmetric_confidence_thresholds) {
    nRatings <- length(grep(pattern = "^theta[0-9]", names(paramDf)))+1
  } else {
    nRatings <- length(grep(pattern = "^thetaUpper[0-9]", names(paramDf)))+1
  }
  vary_sv <-   length(grep(pattern = "^sv[0-9]", names(paramDf), value = T))>1

  if (nConds > 0 ) {
    V <- c(t(paramDf[,paste("v",1:(nConds), sep = "")]))
  } else {
    V <- paramDf$v
    nConds <- 1
  }
  if (vary_sv){
    SV <- c(t((paramDf[,paste("sv",1:(nConds), sep = "")])))
  } else {
    SV <- rep(paramDf$sv, nConds)
  }

  ## Recover confidence thresholds
  if (symmetric_confidence_thresholds) {
    thetas_upper <- c(-1e+32, t(paramDf[,paste("theta",1:(nRatings-1), sep = "")]), 1e+32)
    thetas_lower <- c(-1e+32, t(paramDf[,paste("theta",1:(nRatings-1), sep = "")]), 1e+32)
  } else {
    thetas_upper <- c(-1e+32, t(paramDf[,paste("thetaUpper",1:(nRatings-1), sep = "")]), 1e+32)
    thetas_lower <- c(-1e+32, t(paramDf[,paste("thetaLower",1:(nRatings-1), sep="")]), 1e+32)
  }
  if (thetas_lower[2]>thetas_lower[3]) {
    # For 2DSD the parametrization for lower thetas is different (different confidence scale)
    thetas_lower <- c(-1e+32, rev(thetas_lower[2:(nRatings)]), 1e+32)
    if (symmetric_confidence_thresholds) {
      thetas_lower <- paramDf$a- rev(thetas_upper)
    }
  }
  if (.progress) {
    pb <- progress_bar$new(total = nConds*nRatings*4)
  }
  res <- expand.grid(condition = 1:nConds, stimulus=c("upper", "lower"),
                     response=c("upper","lower"), rating = 1:nRatings) %>%
    group_by(.data$condition, .data$stimulus, .data$response, .data$rating) %>%
    summarise(preddist(.data, thetas_lower = thetas_lower, thetas_upper = thetas_upper,
                paramDf=paramDf, V=V, SV=SV, model=model,
                precision=precision, simult_conf = simult_conf,
                maxrt = maxrt, subdivisions=subdivisions,
                stop.on.error = stop.on.error,
                .progress=.progress, pb=pb))%>%
    mutate(response = if_else(.data$response=="upper", 1, -1),
           stimulus = if_else(.data$stimulus=="upper", 1, -1),
           correct=as.numeric(.data$stimulus==.data$response)) %>%
    ungroup() %>%
    select(c("condition", "stimulus", "response", "correct", "rating", "p", "info", "err"))
  # the last line is to sort the output columns
  # (to combine outputs from predictWEV_Conf and predictRM_Conf)
  res
}



### Predict RT-distribution
#' @rdname predictWEV
#' @export
predictWEV_RT <- function(paramDf, model="dynWEV", precision=1e-5,
                                 maxrt=9, subdivisions = 100L, minrt=NULL,
                                  simult_conf = FALSE,
                                 scaled = FALSE, DistConf=NULL,
                                .progress = TRUE) {
  if (scaled && is.null(DistConf)) {
    message(paste("scaled is TRUE and DistConf is NULL. The rating distribution will",
    " be computed, which will take additional time.", sep=""))
  }
  if (model =="WEVmu") model <- "dynWEV"
  nConds <- length(grep(pattern = "^v[0-9]", names(paramDf), value = T))
  symmetric_confidence_thresholds <- length(grep(pattern = "thetaUpper", names(paramDf), value = T))<1
  if (symmetric_confidence_thresholds) {
    nRatings <- length(grep(pattern = "^theta[0-9]", names(paramDf)))+1
  } else {
    nRatings <- length(grep(pattern = "^thetaUpper[0-9]", names(paramDf)))+1
  }
  vary_sv <-   length(grep(pattern = "^sv[0-9]", names(paramDf), value = T))>1

  if (nConds > 0 ) {
    V <- c(t(paramDf[,paste("v",1:(nConds), sep = "")]))
  } else {
    V <- paramDf$v
    nConds <- 1
  }
  if (vary_sv){
    SV <- c(t((paramDf[,paste("sv",1:(nConds), sep = "")])))
  } else {
    SV <- rep(paramDf$sv, nConds)
  }

  ## Recover confidence thresholds
  if (symmetric_confidence_thresholds) {
    thetas_upper <- c(-1e+32, t(paramDf[,paste("theta",1:(nRatings-1), sep = "")]), 1e+32)
    thetas_lower <- c(-1e+32, t(paramDf[,paste("theta",1:(nRatings-1), sep = "")]), 1e+32)
  } else {
    thetas_upper <- c(-1e+32, t(paramDf[,paste("thetaUpper",1:(nRatings-1), sep = "")]), 1e+32)
    thetas_lower <- c(-1e+32, t(paramDf[,paste("thetaLower",1:(nRatings-1), sep="")]), 1e+32)
  }
  if (thetas_lower[2]>thetas_lower[3]) {
    # For 2DSD the parametrization for lower thetas is different (different confidence scale)
    thetas_lower <- c(-1e+32, rev(thetas_lower[2:(nRatings)]), 1e+32)
    if (symmetric_confidence_thresholds) {
      thetas_lower <- paramDf$a- rev(thetas_upper)
    }
  }

  if (is.null(minrt)) minrt <- paramDf$t0
  df <- expand.grid(rt = seq(minrt, maxrt, length.out = subdivisions),
                    rating = 1:nRatings,
                    response=c("lower", "upper"),
                    stimulus=c(-1, 1),
                    condition = 1:nConds) %>%
    mutate(vth1 = switch(1+as.numeric(model=="2DSD"),
                         if_else(.data$response =="upper", thetas_upper[.data$rating], thetas_lower[(.data$rating)]),
                         if_else(.data$response =="upper", thetas_upper[.data$rating], rev(thetas_lower)[(.data$rating+1)])),
           vth2 = switch(1+as.numeric(model=="2DSD"),
                         if_else(.data$response =="upper", thetas_upper[(.data$rating+1)], thetas_lower[(.data$rating+1)]),
                         if_else(.data$response =="upper", thetas_upper[(.data$rating+1)], rev(thetas_lower)[(.data$rating)])))
  if (model=="dynWEV") {
    dens <- function(df) {
      res <- with(paramDf, dWEV(df$rt, df$vth1,df$vth2,
                    response=as.character(df$response), tau=tau, a=a,
                    v = df$stimulus*V[df$condition],
                    t0 = t0, z = z, sz = sz, st0=st0,
                    sv = SV[df$condition], w=w, svis=svis, sigvis=sigvis,
                    simult_conf = simult_conf, z_absolute = FALSE, precision = precision))
      if (.progress) pb$tick()
      return(data.frame(rt=df$rt, dens=res))
    }
  } else if (model=="2DSD") {
    dens <- function(df) {
      res <- with(paramDf, d2DSD(df$rt, df$vth1,df$vth2,
                   response=as.character(df$response), tau=tau, a=a,
                   v = df$stimulus*V[df$condition],
                   t0 = t0, z = z, sz = sz, st0=st0,
                   sv = SV[df$condition],
                   simult_conf = simult_conf, z_absolute = FALSE, precision = precision))
      if (.progress) pb$tick()
      return(data.frame(rt=df$rt, dens=res))
    }
  }
  if (.progress) {
    pb <- progress_bar$new(total = nConds*nRatings*4)
  }

  df <- df %>% group_by(df[,c("rating", "response", "stimulus", "condition")]) %>%
    summarise(dens(.data))%>%
    mutate(response = if_else(.data$response=="upper", 1, -1))


  if (scaled) {
    ## Scale RT-density to integrate to 1 (for plotting together with simulations)
    # Therefore, divide the density by the probability of a
    # decision-rating-response (as in data.frame DistConf)
    if (is.null(DistConf)) {
      DistConf <- predictWEV_Conf(paramDf, model, precision,
                                  maxrt = maxrt, subdivisions=subdivisions,
                                  simult_conf = simult_conf,
                                  .progress = FALSE) %>%
        ungroup()
    }
    DistConf <- DistConf %>%
      ungroup() %>%
      select(c("rating", "response", "stimulus", "condition", "p"))
    if (is.character(DistConf$response)) {
      DistConf$response <- if_else(DistConf$response=="upper", 1, -1)
    }
    df <- df %>%
      left_join(DistConf, by=c("response", "stimulus", "condition","rating")) %>%
      mutate(densscaled = if_else(.data$p!=0, .data$dens/.data$p, 0)) %>%
      select(-c("p"))
  }
  df <- df %>% mutate(correct = as.numeric(.data$stimulus==.data$response)) %>%
    ungroup() %>%
    select(c("condition", "stimulus", "response", "correct", "rating",
             "rt", "dens", rep("densscaled", as.numeric(scaled))))
  # the last line is to sort the output columns
  # (to combine outputs from predictWEV_RT and predictRM_RT)
  return(df)
}
