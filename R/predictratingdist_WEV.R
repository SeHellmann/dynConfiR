#' Prediction of Confidence Rating and Response Time Distribution in dynaViTE,
#' dynWEV, and 2DSD confidence models
#'
#' \code{predictWEV_Conf} predicts the categorical response distribution of
#' decision and confidence ratings, \code{predictWEV_RT} computes the predicted
#' RT distribution (density) in the 2DSD Model (Pleskac & Busemeyer, 2010) and the
#' dynWEV model (Hellmann et al., 2023), given specific parameter constellations.
#' See \code{\link{ddynaViTE}} and \code{\link{d2DSD}} for more information about parameters.
#'
#' @param paramDf a list or dataframe with one row. Column names should match the names
#' of \link{dynaViTE} and \link{2DSD} model specific parameter names.
#' For different stimulus quality/mean drift rates, names should be `v1`, `v2`, `v3`,....
#' Different `sv` and/or `s` parameters are possible with `sv1`, `sv2`, `sv3`... (`s1`, `s2`, `s3`,...
#' respectively) with equally many steps as for drift rates. Additionally, the confidence
#' thresholds should be given by names with `thetaUpper1`, `thetaUpper2`,..., `thetaLower1`,... or,
#' for symmetric thresholds only by `theta1`, `theta2`,....
#' @param model character scalar. One of "dynaViTE", "dynWEV", or "2DSD".
#' @param precision numerical scalar value. Precision of calculation. Corresponds to the
#' step size of integration w.r.t. `z` and `t0`. Default is 1e-5.
#' @param maxrt numeric. The maximum RT for the integration/density computation.
#' Default: 15 (for `predictWEV_Conf` (integration)), 9 (for `predictWEV_RT`).
#' @param subdivisions integer (default: 100).
#' For \code{predictWEV_Conf} it is used as argument for the inner integral routine.
#' For \code{predictWEV_RT} it is the number of points for which the density is computed.
#' @param minrt numeric or NULL(default). The minimum rt for the density computation.
#' @param  simult_conf logical. Whether in the experiment confidence was reported simultaneously
#' with the decision, as then decision and confidence judgment are assumed to have happened
#' subsequent before response and computations are different, when there is an observable
#' interjudgment time (then `simult_conf` should be FALSE).
#' @param scaled logical. For \code{predictWEV_RT}. Whether the computed density
#' should be scaled to integrate to one (additional column `densscaled`). Otherwise the output
#' contains only the defective density (i.e. its integral is equal to the probability of a
#' response and not 1). If `TRUE`, the argument `DistConf` should be given, if available.
#' Default: `FALSE`.
#' @param DistConf `NULL` or `data.frame`. A `data.frame` or `matrix` with column
#' names, giving the distribution of response and rating choices for
#' different conditions and stimulus categories in the form of the output of
#' \code{predictWEV_Conf}. It is only necessary, if `scaled=TRUE`, because these
#' probabilities are used for scaling. If `scaled=TRUE` and `DistConf=NULL`, it will be computed
#' with the function \code{predictWEV_Conf}, which takes some time and the function will
#' throw a message. Default: `NULL`
#' @param stop.on.error logical. Argument directly passed on to integrate. Default is FALSE,
#' since the densities invoked may lead to slow convergence of the integrals (which are still
#' quite accurate) which causes R to throw an error.
#' @param .progress logical. if TRUE (default) a progress bar is drawn to the console.
#'
#' @return \code{predictWEV_Conf} returns a `data.frame`/`tibble` with columns: `condition`, `stimulus`,
#' `response`, `rating`, `correct`, `p`, `info`, `err`. `p` is the predicted probability of a response
#' and `rating`, given the stimulus category and condition. `info` and `err` refer to the
#' respective outputs of the integration routine used for the computation.
#' \code{predictWEV_RT} returns a `data.frame`/`tibble` with columns: `condition`, `stimulus`,
#' `response`, `rating`, `correct`, `rt` and `dens` (and `densscaled`, if `scaled=TRUE`).
#'
#'
#' @details The function \code{predictWEV_Conf} consists merely of an integration of
#' the response time density, \code{\link{ddynaViTE}} and \code{\link{d2DSD}}, over the response time in a reasonable
#' interval (\code{t0} to `maxrt`). The function \code{predictWEV_RT} wraps these density
#' functions to a parameter set input and a data.frame output.
#' For the argument \code{paramDf}, the output of the fitting function \code{\link{fitRTConf}}
#' with the respective model may be used.
#'
#' @note Different parameters for different conditions are only allowed for drift rate
#' \code{v}, drift rate variability \code{sv}, and process variability `s`. Otherwise, `s` is
#' not required in `paramDf` but set to 1 by default. All other parameters are used for all
#' conditions.
#'
#' @references  Hellmann, S., Zehetleitner, M., & Rausch, M. (2023). Simultaneous modeling of choice, confidence and response time in visual perception. \emph{Psychological Review} 2023 Mar 13. doi: 10.1037/rev0000411. Epub ahead of print. PMID: 36913292.
#'
#' Pleskac, T. J., & Busemeyer, J. R. (2010). Two-Stage Dynamic Signal Detection:
#' A Theory of Choice, Decision Time, and Confidence, \emph{Psychological Review}, 117(3),
#' 864-901. doi:10.1037/a0019737
#'
#'
#' @author Sebastian Hellmann.
#'
#' @name predictWEV
#' @importFrom stats integrate
#' @importFrom progress progress_bar
# @importFrom pracma integral
#' @aliases predictWEV predictdynWEV predict2DSD
#'
#' @examples
#' # Examples for "dynWEV" model (equivalent applicable for "2DSD" model (with less parameters))
#' # 1. Define some parameter set in a data.frame
#' paramDf <- data.frame(a=2.5,v1=0.5, v2=1, t0=0.1,z=0.55,
#'                       sz=0,sv=0.2, st0=0,  tau=3, w=0.3,
#'                       theta1=0.8, svis=0.5, sigvis=0.8)
#'
#' # 2. Predict discrete Choice x Confidence distribution:
#' preds_Conf <- predictWEV_Conf(paramDf, "dynWEV", maxrt = 15)
#' head(preds_Conf)
#' \donttest{
#'   # To set simult_conf=TRUE makes a minor difference in the discrete distribution,
#'   # because we integrate over response times (we just adapt maxrt for comparison)
#'   preds_Conf2 <- predictWEV_Conf(paramDf, "dynWEV", simult_conf = TRUE, maxrt = 15+paramDf$tau)
#'   summary(preds_Conf$p-preds_Conf2$p) # difference in predicted probabilities
#' }
#'
#' # 3. Compute RT density
#' preds_RT <- predictWEV_RT(paramDf, "dynWEV", maxrt=4, subdivisions=200) #(scaled=FALSE)
#' # same output with scaled density column:
#' preds_RT <- predictWEV_RT(paramDf, "dynWEV", maxrt=4, subdivisions=200,
#'                          scaled=TRUE, DistConf = preds_Conf)
#' head(preds_RT)
#' \donttest{
#'   # produces a warning, if scaled=TRUE and DistConf missing
#'   preds_RT <- predictWEV_RT(paramDf, "dynWEV", maxrt=4, subdivisions=200,
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
#' # Use PDFtoQuantiles to get predicted RT quantiles
#' head(PDFtoQuantiles(preds_RT, scaled = FALSE))
#'


#' @rdname predictWEV
#' @export
predictWEV_Conf <- function(paramDf, model="dynaViTE",
                            maxrt=Inf, subdivisions = 100L, simult_conf = FALSE,
                            stop.on.error = FALSE,
                            precision=3,
                            .progress=TRUE){
  if (is.null(model)) {
    if (!("model" %in% names(paramDf))) stop("Either supply model argument or model entry in paramDf argument.")
    model <- paramDf$model
  }
  if (!("lambda" %in% names(paramDf))) {
    if (grepl("dynaViTE|2DSDT", model)) warning("No lambda specified in paramDf. lambda=0 used")
    paramDf$lambda <- 0
  }
  if (grepl("WEVmu|dynWEV|dynaViTE", model)) model <- "dynaViTE"
  if (grepl("2DSD", model)) model <- "2DSD"
  if ("model" %in% names(paramDf)) paramDf$model <- NULL

  paramDf <- fill_optional_params(paramDf, c(st0 = 0, sz = 0, sv = 0))

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
  vary_sv <-   length(grep(pattern = "^sv[0-9]", names(paramDf), value = T))>1
  if (vary_sv){
    SV <- c(t((paramDf[,paste("sv",1:(nConds), sep = "")])))
  } else {
    SV <- rep(paramDf$sv, nConds)
  }
  ## Recover confidence thresholds
  if (symmetric_confidence_thresholds) {
    theta_cols <- grep("^theta[0-9]", names(paramDf), value = TRUE)
    theta_cols <- theta_cols[order(as.integer(sub("^theta", "", theta_cols)))]
    if (length(theta_cols) == 0L) {
      theta_vals <- numeric(0)
    } else {
      theta_vals <- t(paramDf[, theta_cols, drop = FALSE])
    }
    thetas_upper <- c(-1e+32, theta_vals, 1e+32)
    thetas_lower <- c(-1e+32, theta_vals, 1e+32)
  } else {
    theta_upper_cols <- grep("^thetaUpper[0-9]", names(paramDf), value = TRUE)
    theta_upper_cols <- theta_upper_cols[order(as.integer(sub("^thetaUpper", "", theta_upper_cols)))]
    theta_lower_cols <- grep("^thetaLower[0-9]", names(paramDf), value = TRUE)
    theta_lower_cols <- theta_lower_cols[order(as.integer(sub("^thetaLower", "", theta_lower_cols)))]
    if (length(theta_upper_cols) == 0L) {
      theta_upper_vals <- numeric(0)
    } else {
      theta_upper_vals <- t(paramDf[, theta_upper_cols, drop = FALSE])
    }
    if (length(theta_lower_cols) == 0L) {
      theta_lower_vals <- numeric(0)
    } else {
      theta_lower_vals <- t(paramDf[, theta_lower_cols, drop = FALSE])
    }
    thetas_upper <- c(-1e+32, theta_upper_vals, 1e+32)
    thetas_lower <- c(-1e+32, theta_lower_vals, 1e+32)
  }
  if (length(thetas_lower) >= 3 && thetas_lower[2]>thetas_lower[3]) {
    # For 2DSD the parametrization for lower thetas is different (different confidence scale)
    thetas_lower <- c(-1e+32, rev(thetas_lower[2:(nRatings)]), 1e+32)
    if (symmetric_confidence_thresholds) {
      thetas_lower <- paramDf$a- rev(thetas_upper)
    }
  }
  # Because we integrate over the response time, st0 does not matter
  # So, to speed up computations for high values of st0, we set it to 0
  # but add the constant to maxrt
  maxrt <- maxrt + paramDf$st0

  res <- expand.grid(condition = 1:nConds, stimulus=c("upper", "lower"),
                     response=c("upper","lower"), rating = 1:nRatings,
                     p=NA, info=NA, err=NA)
  if (.progress) {
    pb <- progress_bar$new(total = nConds*nRatings*4)
  }


  for (i in 1:nrow(res)) {
    row <- res[i,]
    vth1 <- ifelse(row$response =="upper", thetas_upper[row$rating], thetas_lower[(row$rating)])
    vth2 <- ifelse(row$response =="upper", thetas_upper[(row$rating+1)], thetas_lower[(row$rating+1)])
     if (model == "dynaViTE") {
    p <- integrate(function(rt) return(ddynaViTE(rt, response=as.character(row$response),
                                      vth1,vth2, a=paramDf$a,
                                      v = (-1)^(row$stimulus=="lower")*V[row$condition],
                                      t0 = paramDf$t0, z = paramDf$z, sz = paramDf$sz, sv = SV[row$condition],
                                      st0=0, # we integrate over t, so this does not change results,
                                      # but speeds up computations considerably
                                      tau=paramDf$tau,
                                       w=paramDf$w, svis=paramDf$svis, sigvis=paramDf$sigvis,
                                      lambda=paramDf$lambda, s = S[row$condition],
                                      simult_conf = simult_conf,
                                      z_absolute = FALSE, precision = precision)),
              lower=paramDf$t0, upper=maxrt, subdivisions = subdivisions,
                 stop.on.error = stop.on.error)
  } else if (model=="2DSD") {
    p <- integrate(function(rt) return(d2DSD(rt, response=as.character(row$response),
                                             vth1,vth2, a=paramDf$a,
                                             v = (-1)^(row$stimulus=="lower")*V[row$condition],
                                             t0 = paramDf$t0, z = paramDf$z, sz = paramDf$sz,
                                             sv = SV[row$condition], st0=0, tau=paramDf$tau,
                                             lambda=paramDf$lambda, s = S[row$condition],
                                             simult_conf = simult_conf,
                                             z_absolute = FALSE, precision = precision)),
                   lower=paramDf$t0, upper=maxrt, subdivisions = subdivisions,
                   stop.on.error = stop.on.error)
  } else { stop("model must contain dynaViTE, dynWEV, or 2DSD") }
  if (.progress) pb$tick()
  res[i, 5:7] <- list(p = p$value, info = p$message, err = p$abs.error)
  }

  res$response <- 2*as.numeric(res$response=="upper")-1
  res$stimulus <- 2*as.numeric(res$stimulus=="upper")-1
  res$correct <- as.numeric(res$stimulus==res$response)
  res <- res[c("condition", "stimulus", "response", "correct", "rating", "p", "info", "err")]
  # the last line is to sort the output columns
  # (to combine outputs from predictWEV_Conf and predictRM_Conf)
  res
}



### Predict RT-distribution
#' @rdname predictWEV
#' @export
predictWEV_RT <- function(paramDf, model=NULL,
                                 maxrt=9, subdivisions = 100L, minrt=NULL,
                                  simult_conf = FALSE,
                                 scaled = FALSE, DistConf=NULL, precision=3,
                                .progress = TRUE) {
  if (scaled && is.null(DistConf)) {
    message(paste("scaled is TRUE and DistConf is NULL. The rating distribution will",
    " be computed, which will take additional time.", sep=""))
  }
  if (is.null(model)) {
    if (!("model" %in% names(paramDf))) stop("Either supply model argument or model entry in paramDf argument.")
    model <- paramDf$model
  }
  if (!("lambda" %in% names(paramDf))) {
    if (grepl("dynaViTE|2DSDT", model)) warning("No lambda specified in paramDf. lambda=0 used")
    paramDf$lambda <- 0
  }
  if (grepl("WEVmu|dynWEV|dynaViTE", model)) model <- "dynaViTE"
  if (grepl("2DSD", model)) model <- "2DSD"
  if ("model" %in% names(paramDf)) paramDf$model <- NULL

  paramDf <- fill_optional_params(paramDf, c(st0 = 0, sz = 0, sv = 0))

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
  vary_sv <-   length(grep(pattern = "^sv[0-9]", names(paramDf), value = T))>1
  if (vary_sv){
    SV <- c(t((paramDf[,paste("sv",1:(nConds), sep = "")])))
  } else {
    SV <- rep(paramDf$sv, nConds)
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
    theta_cols <- grep("^theta[0-9]", names(paramDf), value = TRUE)
    theta_cols <- theta_cols[order(as.integer(sub("^theta", "", theta_cols)))]
    if (length(theta_cols) == 0L) {
      theta_vals <- numeric(0)
    } else {
      theta_vals <- t(paramDf[, theta_cols, drop = FALSE])
    }
    thetas_upper <- c(-1e+32, theta_vals, 1e+32)
    thetas_lower <- c(-1e+32, theta_vals, 1e+32)
  } else {
    theta_upper_cols <- grep("^thetaUpper[0-9]", names(paramDf), value = TRUE)
    theta_upper_cols <- theta_upper_cols[order(as.integer(sub("^thetaUpper", "", theta_upper_cols)))]
    theta_lower_cols <- grep("^thetaLower[0-9]", names(paramDf), value = TRUE)
    theta_lower_cols <- theta_lower_cols[order(as.integer(sub("^thetaLower", "", theta_lower_cols)))]
    if (length(theta_upper_cols) == 0L) {
      theta_upper_vals <- numeric(0)
    } else {
      theta_upper_vals <- t(paramDf[, theta_upper_cols, drop = FALSE])
    }
    if (length(theta_lower_cols) == 0L) {
      theta_lower_vals <- numeric(0)
    } else {
      theta_lower_vals <- t(paramDf[, theta_lower_cols, drop = FALSE])
    }
    thetas_upper <- c(-1e+32, theta_upper_vals, 1e+32)
    thetas_lower <- c(-1e+32, theta_lower_vals, 1e+32)
  }
  if (length(thetas_lower) >= 3 && thetas_lower[2]>thetas_lower[3]) {
    # For 2DSD the parametrization for lower thetas is different (different confidence scale)
    thetas_lower <- c(-1e+32, rev(thetas_lower[2:(nRatings)]), 1e+32)
    if (symmetric_confidence_thresholds) {
      thetas_lower <- paramDf$a- rev(thetas_upper)
    }
  }
  if (!("lambda" %in% names(paramDf))) paramDf$lambda <- 0

  if (is.null(minrt)) minrt <- paramDf$t0
  rt = seq(minrt, maxrt, length.out = subdivisions)
  df <- expand.grid(rt = rt,
                    rating = 1:nRatings,
                    response=c("upper", "lower"),
                    stimulus=c("upper", "lower"),
                    condition = 1:nConds, dens=NA)
  if (scaled) {
    ## Scale RT-density to integrate to 1 (for plotting together with simulations)
    # Therefore, divide the density by the probability of a
    # decision-rating-response (as in data.frame DistConf)
    if (is.null(DistConf)) {
      DistConf <- predictWEV_Conf(paramDf, model, precision,
                                  maxrt = maxrt, subdivisions=subdivisions,
                                  simult_conf = simult_conf,
                                  .progress = FALSE)
    }
    DistConf <- DistConf[,c("rating", "response", "stimulus", "condition", "p")]
    if (!is.character(DistConf$response)) {
      DistConf[DistConf$response==1, "response"] <- "upper"
      DistConf[DistConf$response==-1, "response"] <- "lower"
    }
    if (!is.character(DistConf$stimulus)) {
      DistConf[DistConf$stimulus==1, "stimulus"] <- "upper"
      DistConf[DistConf$stimulus==-1, "stimulus"] <- "lower"
    }
    df$densscaled <- NA
  }


  if (.progress) {
    pb <- progress_bar$new(total = nConds*nRatings*4)
  }
  for ( i in 1:(nRatings*2*2*nConds)) {
    cur_row <- df[1+((i-1)*subdivisions),]
    vth1 = ifelse(cur_row$response =="upper", thetas_upper[cur_row$rating], thetas_lower[(cur_row$rating)])
    vth2 = ifelse(cur_row$response =="upper", thetas_upper[(cur_row$rating+1)], thetas_lower[(cur_row$rating+1)])
    v <- V[cur_row$condition]*(-1)^(cur_row$stimulus=="lower")


    if (model=="dynaViTE") {
      df[(1:subdivisions) + subdivisions*(i-1), "dens"] <-
        ddynaViTE(rt, response=as.character(cur_row$response),
             vth1,vth2,v = v,
             tau=paramDf$tau, a=paramDf$a,
             t0 = paramDf$t0, z = paramDf$z, sz = paramDf$sz, st0=paramDf$st0,
             sv = SV[cur_row$condition], lambda=paramDf$lambda,s = S[cur_row$condition],
             w=paramDf$w, svis=paramDf$svis, sigvis=paramDf$sigvis,
             simult_conf = simult_conf, z_absolute = FALSE, precision = precision)
    } else if (model=="2DSD") {
      df[(1:subdivisions) + subdivisions*(i-1), "dens"] <-
        d2DSD(rt, response=as.character(cur_row$response),
              vth1,vth2,v = v,
              tau=paramDf$tau, a=paramDf$a,
              t0 = paramDf$t0, z = paramDf$z, sz = paramDf$sz, st0=paramDf$st0,
              sv = SV[cur_row$condition], lambda=paramDf$lambda,s = S[cur_row$condition],
              simult_conf = simult_conf, z_absolute = FALSE, precision = precision)
    } else { stop("model must contain dynaViTE, dynWEV, or 2DSD") }

    if (scaled) {
      P <- DistConf[DistConf$condition==cur_row$condition &
                      DistConf$response==cur_row$response &
                      DistConf$rating == cur_row$rating &
                      DistConf$stimulus==cur_row$stimulus,]$p
      if (P != 0) {
        df[(1:subdivisions) + subdivisions*(i-1), "densscaled"] <-
          df[(1:subdivisions) + subdivisions*(i - 1), "dens"]/P
      } else {
        df[(1:subdivisions) + subdivisions*(i-1), "densscaled"] <- 0
      }
    }
    if (.progress) pb$tick()
  }
  df$response <- 2*as.numeric(df$response=="upper")-1
  df$stimulus <- 2*as.numeric(df$stimulus=="upper")-1
  df$correct <-  as.numeric(df$stimulus==df$response)
  df <- df[,c("condition", "stimulus", "response", "correct", "rating",
              "rt", "dens", rep("densscaled", as.numeric(scaled)))]
  # the last line is to sort the output columns
  # (to combine outputs from predictWEV_RT and predictDDConf_RT)
  return(df)
}
