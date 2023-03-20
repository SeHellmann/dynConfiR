#' Prediction of Confidence Rating and Reaction Time Distribution in race models of confidence
#'
#' \code{predictRM_Conf} predicts the categorical response distribution of
#' decision and confidence ratings, \code{predictRM_RT} computes the
#' RT distribution (density) in the independent and partially anti-correlated
#' race models  (Hellmann et al., 2023), given specific parameter
#' constellations. See \link{RaceModels} for more information about the models
#' and parameters.
#'
#' @param paramDf a list or data frame with one row. Column names should match the names of
#' \link{RaceModels} parameter names (only `mu1` and `mu2` are not used in this context but
#' replaced by the parameter `v`). For different stimulus quality/mean
#' drift rates, names should be `v1`, `v2`, `v3`,....
#' Different `s` parameters are possible with `s1`, `s2`, `s3`,... with equally many steps as for drift rates. Additionally, the confidence
#' thresholds should be given by names with `thetaUpper1`, `thetaUpper2`,..., `thetaLower1`,... or,
#' for symmetric thresholds only by `theta1`, `theta2`,....
#' @param model character scalar. One of "IRM" or "PCRM". ("IRMt" and "PCRMt" will also be accepted. In that case,
#' time_scaled is set to TRUE.)
#' @param time_scaled logical. Whether the confidence measure should be scaled by 1/sqrt(rt). Default: FALSE.
#' (It is set to TRUE, if model is "IRMt" or "PCRMt")
#' @param maxrt numeric. The maximum RT for the
#' integration/density computation. Default: 15 (for `predictRM_Conf`
#' (integration)), 9 (for `predictRM_RT`).
#' @param subdivisions \code{integer} (default: 100).
#' For \code{predictRM_Conf} it is used as argument for the inner integral routine.
#' For \code{predictRM_RT} it is the number of points for which the density is computed.
#' @param minrt numeric or NULL(default). The minimum rt for the density computation.
#' @param scaled logical. For \code{predictRM_RT}. Whether the computed density
#' should be scaled to integrate to one (additional column `densscaled`). Otherwise the output
#' contains only the defective density (i.e. its integral is equal to the probability of a
#' response and not 1). If `TRUE`, the argument `DistConf` should be given, if available.
#' Default: `FALSE`.
#' @param DistConf `NULL` or `data.frame`. A `data.frame` or `matrix` with column
#' names, giving the distribution of response and rating choices for
#' different conditions and stimulus categories in the form of the output of
#' \code{predictRM_Conf}. It is only necessary, if `scaled=TRUE`, because these
#' probabilities are used for scaling. If `scaled=TRUE` and `DistConf=NULL`, it will be computed
#' with the function \code{predictRM_Conf}, which takes some time and the function will
#' throw a message. Default: `NULL`
#' @param stop.on.error logical. Argument directly passed on to integrate. Default is FALSE,
#' since the densities invoked may lead to slow convergence of the integrals (which are still
#' quite accurate) which causes R to throw an error.
#' @param .progress logical. If TRUE (default) a progress bar is drawn to the console.
#'
#' @return \code{predictRM_Conf} returns a `data.frame`/`tibble` with columns: `condition`, `stimulus`,
#' `response`, `rating`, `correct`, `p`, `info`, `err`. `p` is the predicted probability of a response
#' and `rating`, given the stimulus category and condition. `info` and `err` refer to the
#' respective outputs of the integration routine used for the computation.
#' \code{predictRM_RT} returns a `data.frame`/`tibble` with columns: `condition`, `stimulus`,
#' `response`, `rating`, `correct`, `rt` and `dens` (and `densscaled`, if `scaled=TRUE`).
#'
#'
#' @details The function \code{predictRM_Conf} consists merely of an integration of
#' the response time density, \code{\link{dIRM}} and \code{\link{dPCRM}}, over the
#' response time in a reasonable interval (0 to `maxrt`). The function
#' \code{predictRM_RT} wraps these density
#' functions to a parameter set input and a data.frame output.
#' For the argument \code{paramDf}, the output of the fitting function \code{\link{fitRTConf}}
#' with the respective model may be used.
#'
#' The drift rate parameters differ from those used in \code{\link{dIRM}}/\code{\link{dPCRM}}
#' since in many perceptual decision experiments the drift on one accumulator is assumed to
#' be the negative of the other. The drift rate of the correct accumulator is `v` (`v1`, `v2`,
#' ... respectively) in `paramDf`.
#'
#' @note Different parameters for different conditions are only allowed for drift rate,
#' \code{v}, and process variability \code{s}. All other parameters are used for all
#' conditions.
#'
#' @references  Hellmann, S., Zehetleitner, M., & Rausch, M. (2023). Simultaneous modeling of choice, confidence and response time in visual perception. \emph{Psychological Review} 2023 Mar 13. doi: 10.1037/rev0000411. Epub ahead of print. PMID: 36913292.
#'
#' @author Sebastian Hellmann.
#'
#' @name predictRM
#' @importFrom stats integrate
#' @import dplyr
#' @importFrom progress progress_bar
#' @importFrom magrittr %>%
#' @importFrom rlang .data
# @importFrom pracma integral
#' @aliases predictIRM predictPCRM
#' @importFrom Rcpp evalCpp
#'
#' @examples
#' # Examples for "PCRM" model (equivalent applicable for "IRM" model)
#' # 1. Define some parameter set in a data.frame
#' paramDf <- data.frame(a=2,b=2, v1=0.5, v2=1, t0=0.1,st0=0,
#'                       wx=0.6, wint=0.2, wrt=0.2,
#'                       theta1=4)
#'
#' # 2. Predict discrete Choice x Confidence distribution:
#' preds_Conf <- predictRM_Conf(paramDf, "PCRM", time_scaled=TRUE)
#' # equivalent:
#' preds_Conf <- predictRM_Conf(paramDf, "PCRMt")
#' head(preds_Conf)
#'
#' # 3. Compute RT density
#' preds_RT <- predictRM_RT(paramDf, "PCRMt", maxrt=7, subdivisions=50)
#' # same output with scaled density column:
#' preds_RT <- predictRM_RT(paramDf, "PCRMt", maxrt=7, subdivisions=50,
#'                          scaled=TRUE, DistConf = preds_Conf)
#' head(preds_RT)
#' \donttest{
#'   # produces a warning, if scaled=TRUE and DistConf missing
#'   preds_RT <- predictRM_RT(paramDf, "PCRMt", maxrt=7, subdivisions=50,
#'                            scaled=TRUE)
#' }
#'
#' \donttest{
#'   # Example of visualization
#'   library(ggplot2)
#'   preds_Conf$rating <- factor(preds_Conf$rating, labels=c("unsure", "sure"))
#'   preds_RT$rating <- factor(preds_RT$rating, labels=c("unsure", "sure"))
#'   ggplot(preds_Conf, aes(x=interaction(rating, response), y=p))+
#'     geom_bar(stat="identity")+
#'     facet_grid(cols=vars(stimulus), rows=vars(condition), labeller = "label_both")
#'   ggplot(preds_RT, aes(x=rt, color=interaction(rating, response), y=dens))+
#'     geom_line(stat="identity")+
#'     facet_grid(cols=vars(stimulus), rows=vars(condition), labeller = "label_both")+
#'     theme(legend.position = "bottom")
#'   ggplot(aggregate(densscaled~rt+correct+rating+condition, preds_RT, mean),
#'          aes(x=rt, color=rating, y=densscaled))+
#'     geom_line(stat="identity")+
#'     facet_grid(cols=vars(condition), rows=vars(correct), labeller = "label_both")+
#'     theme(legend.position = "bottom")
#' }
#' \donttest{
#'   # Use PDFtoQuantiles to get predicted RT quantiles
#'   # (produces warning because of few rt steps (--> inaccurate calculations))
#'   PDFtoQuantiles(preds_RT, scaled = FALSE)
#' }
#'


#' @rdname predictRM
#' @export
predictRM_Conf <- function(paramDf, model="IRM", time_scaled = FALSE,
                                maxrt=15, subdivisions = 100L, stop.on.error=FALSE,
                                .progress=TRUE){
  #### Check model argument
  if (model=="IRMt") {
    time_scaled=TRUE
  }
  if (model=="PCRMt") {
    time_scaled=TRUE
  }
  if (!model %in% c("IRM", "PCRM","IRMt", "PCRMt")) stop("model must be 'IRM', 'PCRM', 'IRMt' or 'PCRMt'")

  if (!time_scaled) {
    paramDf$wx = 1
    paramDf$wrt = 0
    paramDf$wint = 0
  }

  nConds <- length(grep(pattern = "^v[0-9]", names(paramDf), value = T))
  symmetric_confidence_thresholds <- length(grep(pattern = "thetaUpper", names(paramDf), value = T))<1
  if (symmetric_confidence_thresholds) {
    nRatings <- length(grep(pattern = "^theta[0-9]", names(paramDf)))+1
  } else {
    nRatings <- length(grep(pattern = "^thetaUpper[0-9]", names(paramDf)))+1
  }

  if (nConds > 0 ) {
    V <- c(t(paramDf[,paste("v",1:(nConds), sep = "")]))
  } else {
    V <- paramDf$v
    nConds <- 1
  }
  vary_s <-   length(grep(pattern = "^s[0-9]", names(paramDf), value = T))>1
  if (vary_s){
    S <- c(t((paramDf[,paste("s",1:(nConds), sep = "")])))
  } else {
    if ("s" %in% names(paramDf)) {
      S <- rep(paramDf$s, nConds)
    } else {
      S <- rep(1, nConds)
    }
  }
  ## Recover confidence thresholds
  if (symmetric_confidence_thresholds) {
    if (nRatings==1) {
      thetas_1 <- c(1e-32, 1e+64)
      thetas_2 <- c(1e-32, 1e+64)
    } else {
      thetas_1 <- c(1e-32, t(paramDf[,paste("theta",1:(nRatings-1), sep = "")]), 1e+64)
      thetas_2 <- c(1e-32, t(paramDf[,paste("theta",1:(nRatings-1), sep = "")]), 1e+64)
    }
  } else {
    thetas_1 <- c(1e-32, t(paramDf[,paste("thetaLower",1:(nRatings-1), sep = "")]), 1e+64)
    thetas_2 <- c(1e-32, t(paramDf[,paste("thetaUpper",1:(nRatings-1), sep="")]), 1e+64)
  }
  # Because we integrate over the response time, st0 does not matter
  # So, to speed up computations for high values of st0, we set it to 0
  # but add the constant to maxrt
  maxrt <- maxrt + paramDf$st0
  paramDf$st0 <- 0

  if (.progress) {
    pb <- progress_bar$new(total = nConds*nRatings*4)
  }

  if (model=="IRM") {
    help_fct <- function(row) {
      p <- integrate(function(rt) return(dIRM(rt, response=row$response,
                                              mu1=row$mu1, mu2 = row$mu2, s=row$s,
                                              a = row$a, b=row$b,
                                              wx = row$wx, wrt= row$wrt, wint = row$wint,
                                              th1 = row$th1, th2=row$th2,
                                              t0  = row$t0,
                                              st0 = row$st0, time_scaled=time_scaled)),
                     lower=paramDf$t0, upper=maxrt, subdivisions = subdivisions,
                     stop.on.error = stop.on.error)
      if (.progress) pb$tick()
      return(data.frame(p = p$value, info = p$message, err = p$abs.error))
    }
  } else {
    help_fct <- function(row) {
      p <- integrate(function(rt) return(dPCRM(rt, response=row$response,
                                               mu1=row$mu1, mu2 = row$mu2, s=row$s,
                                               a = row$a, b=row$b,
                                               wx = row$wx, wrt= row$wrt, wint = row$wint,
                                               th1 = row$th1, th2=row$th2,
                                               t0  = row$t0,
                                               st0 = row$st0, time_scaled=time_scaled)),
                     lower=paramDf$t0, upper=maxrt, subdivisions = subdivisions,
                     stop.on.error = stop.on.error)
      if (.progress) pb$tick()
      return(data.frame(p = p$value, info = p$message, err = p$abs.error))
    }
  }


  res <- expand.grid(condition = 1:nConds, stimulus=c(1,2),
                     response=c(1,2), rating = 1:nRatings) %>%
    mutate(a    = paramDf$a,
           b    = paramDf$b,
           s    = S[.data$condition],
           th1  = if_else(.data$response == 1, thetas_1[(.data$rating)], thetas_2[(.data$rating)]),
           th2  = if_else(.data$response == 1, thetas_1[(.data$rating+1)], thetas_2[(.data$rating + 1)]),
           mu1  = V[.data$condition]*(-1)^(1+.data$stimulus),
           mu2  = V[.data$condition]*(-1)^(.data$stimulus),
           wx   = paramDf$wx,
           wrt  = paramDf$wrt,
           wint = paramDf$wint,
           t0   = paramDf$t0,
           st0  = paramDf$st0) %>%
    group_by(.data$condition, .data$stimulus, .data$response, .data$rating) %>%
    summarise(help_fct(cur_data_all()))%>%
    mutate(correct=as.numeric(.data$stimulus==.data$response)) %>%
    ungroup() %>%
    select(c("condition", "stimulus", "response", "correct", "rating", "p", "info", "err"))
  # the last line is to sort the output columns
  # (to combine outputs from predictWEV_Conf and predictRM_Conf)
  res
}



### Predict RT-distribution
#' @rdname predictRM
#' @export
predictRM_RT <- function(paramDf, model="IRM", time_scaled = FALSE,
                                 maxrt=9, subdivisions = 100L,  minrt=NULL,
                                 scaled = FALSE, DistConf=NULL,
                                .progress = TRUE) {
  #### Check model argument
  if (model=="IRMt") {
    time_scaled=TRUE
  }
  if (model=="PCRMt") {
    time_scaled=TRUE
  }
  if (!time_scaled) {
    paramDf$wx = 1
    paramDf$wrt = 0
    paramDf$wint = 0
  }

  if (scaled && is.null(DistConf)) {
    message(paste("scaled is TRUE and DistConf is NULL. The rating distribution will",
    " be computed, which will take additional time.", sep=""))
  }
  nConds <- length(grep(pattern = "^v[0-9]", names(paramDf), value = T))
  symmetric_confidence_thresholds <- length(grep(pattern = "thetaUpper", names(paramDf), value = T))<1
  if (symmetric_confidence_thresholds) {
    nRatings <- length(grep(pattern = "^theta[0-9]", names(paramDf)))+1
  } else {
    nRatings <- length(grep(pattern = "^thetaUpper[0-9]", names(paramDf)))+1
  }

  if (nConds > 0 ) {
    V <- c(t(paramDf[,paste("v",1:(nConds), sep = "")]))
  } else {
    V <- paramDf$v
    nConds <- 1
  }
  vary_s <-   length(grep(pattern = "^s[0-9]", names(paramDf), value = T))>1
  if (vary_s){
    S <- c(t((paramDf[,paste("s",1:(nConds), sep = "")])))
  } else {
    if ("s" %in% names(paramDf)) {
      S <- rep(paramDf$s, nConds)
    } else {
      S <- rep(1, nConds)
    }
  }


  ## Recover confidence thresholds
  if (symmetric_confidence_thresholds) {
    if (nRatings==1) {
      thetas_1 <- c(1e-32, 1e+32)
      thetas_2 <- c(1e-32, 1e+32)
    } else {
      thetas_1 <- c(1e-32, t(paramDf[,paste("theta",1:(nRatings-1), sep = "")]), 1e+64)
      thetas_2 <- c(1e-32, t(paramDf[,paste("theta",1:(nRatings-1), sep = "")]), 1e+64)
    }
  } else {
    thetas_1 <- c(1e-32, t(paramDf[,paste("thetaLower",1:(nRatings-1), sep = "")]), 1e+64)
    thetas_2 <- c(1e-32, t(paramDf[,paste("thetaUpper",1:(nRatings-1), sep="")]), 1e+64)
  }

  if (is.null(minrt)) minrt <- paramDf$t0
  df <- expand.grid(rt = seq(minrt, maxrt, length.out = subdivisions),
                    rating = 1:nRatings,
                    response=c(1,2),
                    stimulus=c(1,2),
                    condition = 1:nConds) %>%
    mutate(th1 = if_else(.data$response==1, thetas_1[(.data$rating)], thetas_2[(.data$rating)]),
           th2 = if_else(.data$response==1, thetas_1[(.data$rating+1)], thetas_2[(.data$rating+1)]),
           mu1 = V[.data$condition]*(-1)^(1+.data$stimulus),
           mu2 = V[.data$condition]*(-1)^(.data$stimulus),
           s = S[.data$condition])
  if (.progress) {
    pb <- progress_bar$new(total = nConds*nRatings*4)
  }
  if (grepl("IRM", model)) {
    dens <- function(df) {
      res <- dIRM(df$rt, df$response[1],
                  mu1=df$mu1[1], mu2 = df$mu2[1], s=df$s[1],
                  a = paramDf$a, b=paramDf$b,
                  wx = paramDf$wx, wrt= paramDf$wrt, wint = paramDf$wint,
                  th1 = df$th1[1], th2=df$th2[1],
                  t0  = paramDf$t0,
                  st0 = paramDf$st0, time_scaled=time_scaled)
      if (.progress) pb$tick()
      return(data.frame(rt=df$rt, dens=res))
    }
  } else if (grepl("PCRM", model)) {
    dens <- function(df) {
      res <- dPCRM(df$rt, df$response[1],
                   mu1=df$mu1[1], mu2 = df$mu2[1], s=df$s[1],
                   a = paramDf$a, b=paramDf$b,
                   wx = paramDf$wx, wrt= paramDf$wrt, wint = paramDf$wint,
                   th1 = df$th1[1], th2=df$th2[1],
                   t0  = paramDf$t0,
                   st0 = paramDf$st0, time_scaled=time_scaled)
      if (.progress) pb$tick()
      return(data.frame(rt=df$rt, dens=res))
    }
  } else { stop("model must be IRM, IRMt, PCRM or PCRMt") }

  df <- df %>% group_by(df[,c("rating", "response", "stimulus", "condition")]) %>%
    summarise(dens(.data))


  if (scaled) {
    ## Scale RT-density to integrate to 1 (for plotting together with simulations)
    # Therefore, divide the density by the probability of a
    # decision-rating-response (as in data.frame DistConf)
    if (is.null(DistConf)) {
      DistConf <- predictRM_Conf(paramDf, model, time_scaled,
                                 maxrt = maxrt, subdivisions=subdivisions,
                                 .progress = FALSE) %>%
        ungroup()
    }
    DistConf <- DistConf %>%
      ungroup() %>%
      select(c("rating", "response", "stimulus", "condition", "p"))
    if (is.character(DistConf$response)) {
      DistConf$response <- as.integer(as.factor(DistConf$response))
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
