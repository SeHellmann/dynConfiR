#' Simulation of confidence ratings and RTs in dWEV and 2DSD confidence models
#'
#' Simulates the decision responses and reaction times together with a discrete confidence judgement  in the 2DSD Model (Pleskac & Busemeyer, 2010) and the dWEV model (Hellmann & Rausch), given specific parameter constellations.
#' See \code{\link{dWEV}} and \code{\link{d2DSD}} for more information about parameters. Also computes the Gamma rank correlation between the confidence ratings and condition (task difficulty), reaction times and accuracy in the simulated output.
#'
#' @param n integer. The number of samples (per condition and stimulus direction) generated. Total number of samples is \code{n*nConditions*length(stimulus)}.
#' @param paramDf a list or dataframe with one row. Column names should match the names of dWEV and 2DSD model specific parameter names. For different stimulus quality/mean drift rates, names should be v1, v2, v3,.... Different sv parameters are possible with sv1, sv2, sv3... with equally many steps as for drift rates. Additionally, the confidence thresholds should be given by names with thetaUpper1, thetaUpper2,..., thetaLower1,... or, for symmetric thresholds only by theta1, theta2,....
#' @param model character scalar. One of "WEVmu", "dynWEV", or "2DSD".
#' @param delta numeric. Discretization steps for simulations with the stochastic process (used, if method!="rtdists")
#' @param maxrt numeric. Maximum reaction time returned, if method!="rtdists". If the simulation of the stochastic process exceeds a rt of maxrt,
#' the response will be set to 0 and the maxrt will be returned as rt.
#' @param  simult_conf logical. TRUE, if in the experiment confidence was reported simultaneously
#' with the decision, as then decision and confidence judgment are assumed to have happened
#' subsequent before response and tau is added to the simulated decision time. If FALSE,
#' returned response time will only be decision time plus non-judgment time component.
#' @param gamma logical.
#' @param agg_simus logical. Simulation is done on a trial basis with rt's outcome. If TRUE, the simulations will be aggregated over RTs to return only the distribution of response and confidence ratings. Default: FALSE.
#'
#' @param stimulus numeric vector. Either 1, -1 or c(-1, 1) (default). Together with condition represents the experimental situation. In a 2AFC task the presented stimulus belongs to one of two categories. In the default setting trials with
#' both categories presented are simulated but one can choose to simulate only trials with the stimulus coming from one category (1 for the category that is associated with positive drift in the decision process where "upper"/1 responses are considered correct and -1 correspondingly for negative drifts and "lower"/-1 correct decisions).
#' @param method character. Method for the simulation. If "rtdists", simulations are done using the function rdiffusion from the rtdists package. Else (default), for each observation the
#' stochastic process is simulated in discrete steps.
#' @param precision \code{numerical} scalar value. Precision of calculation. Corresponds roughly to the number of decimals of the predicted CDFs that are calculated accurately. Default is 3. This argument is given directly to \code{\link[rtdists:Diffusion]{rdiffusion}} used for generating samples in the decision process.
#' @param seed numerical. Seeding for non-random data generation.
#'
#' @return Depending on gamma and agg_simus. If gamma is TRUE, returns a list with elements:
#' "simus" (the simulated data frame) and "gamma", which is again a list with elements
#' "condition", "rt" and "correct", each a tibble with two columns (see details for more
#' information). If gamma is FALSE, returns a data frame with columns: condition, stimulus,
#' response, correct, rt, conf (the continuous confidence measure) and rating (the discrete
#' confidence rating) or (if agg_simus=TRUE): condition, stimulus, response, correct, rating
#' and p (for the probability of a response and rating, given the condition and stimulus).
#'
#'
#' @details The function combines the random generator \code{\link[rtdists:Diffusion]{rdiffusion}} from the package \code{rtdists} for the decision process outputs and \code{rnorm} to produce the confidence measure in the respective model.
#' The confidence outputs are then binned according to the given thresholds. The output of the fitting function \code{\link{fitRTConf}} with the respective model fits the argument paramDf for simulation.
#' The Gamma coefficients are computed seperately for correct/incorrect responses for the correlation of confidence ratings with condition and rt and seperately for conditions for the correlation of accuracy and confidence. The
#' resulting tibbles in the output thus have two columns. One for the grouping variable and one for the Gamma coefficient.
#'
#' @note Different parameters for different conditions are only allowed for drift rate, \code{v}, and drift rate variability, \code{sv}. All other parameters are used for all conditions.
#'
#' @references Pleskac, T. J., & Busemeyer, J. R. (2010). Two-Stage Dynamic Signal Detection: A Theory of Choice, Decision Time, and Confidence, \emph{Psychological Review}, 117(3), 864-901. doi:10.1037/a0019737
#'
#' Rausch, M., Hellmann, S., & Zehetleitner, M. (2018). Confidence in masked orientation judgments is informed by both evidence and visibility. \emph{Attention, Perception, & Psychophysics}, 80(1), 134–154.  doi: 10.3758/s13414-017-1431-5
#'
#'
#' @author Sebastian Hellmann.
#'
#' @name simulateWEV
#' @importFrom rtdists rdiffusion
#' @importFrom stats rnorm
#' @import dplyr
#' @importFrom magrittr %>%
#' @importFrom Hmisc rcorr.cens
#' @importFrom rlang .data
# @importFrom pracma integral
#' @aliases rWEV r2DSD simulate2DSD
#'

## When given vectorised parameters, n is the number of replicates for each parameter set
#' @rdname simulateWEV
#' @export
rWEV <- function (paramDf, n=1e+4,  model = "dynWEV", simult_conf = TRUE, gamma = FALSE, agg_simus=FALSE,
                         stimulus = c(-1,1), method = "Cpp",  precision = 3, delta=0.01, maxrt=15, seed=NULL)
{
  if (!is.null(seed)) {
    set.seed(seed)
  }
  if (model=="WEVmu") model <- "dynWEV"
  if (!(model %in% c("dynWEV", "2DSD"))) stop("Only models dynWEV (alias: WEVmu) and 2DSD are allowed!")

  if (!(all(stimulus %in% c(-1, 1)))) {
    stop(paste("Not accepted value for stimulus: ", paste(stimulus, collapse=", "),". Must be either 1, -1 or c(-1,1).", sep=""))
  }


  ## recover parameters from paramDf
  a <- paramDf$a
  z <- paramDf$z
  sz <- paramDf$sz
  t0 <- paramDf$t0
  st0 <- paramDf$st0
  tau = paramDf$tau

  nConds <- length(grep(pattern = "^v[0-9]", names(paramDf), value = T))
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

  symmetric_confidence_thresholds <- length(grep(pattern = "thetaUpper", names(paramDf), value = T))<1
  if (symmetric_confidence_thresholds) {
    nRatings <- length(grep(pattern = "^theta[0-9]", names(paramDf)))+1
  } else {
    nRatings <- length(grep(pattern = "^thetaUpper[0-9]", names(paramDf)))+1
  }

  if (method == "rtdists") {
    simus <- expand.grid(condition = rep(1:nConds, each=n), stimulus=stimulus)
    ### Produce decision responses and confidence measures according to the given model
    if (model == "2DSD") {
      if (all(SV==0)) {
        simus <- cbind(simus, rdiffusion(n=nrow(simus), a=rep(a, nrow(simus)), v=V[simus$condition]*simus$stimulus, t0 = t0,
                                     z = z*a, d = 0, sz = a*sz, sv = 0,
                                     st0=st0, s=1,
                                     precision = precision))
        simus$conf <- rnorm(n=nrow(simus),
                               mean=a * 0^(simus$response=="lower") +tau*V[simus$condition]*simus$stimulus,
                               sd=sqrt(tau))
      } else {
        simus$d <- rnorm(n=nrow(simus),
                       mean=V[simus$condition]*simus$stimulus,
                       sd = SV[simus$condition])
        simus <- cbind(simus,
                     rdiffusion(n=nrow(simus), a=rep(a, nrow(simus)), v=simus$d, t0 = t0,
                                z = z*a, d = 0, sz = a*sz, sv = 0,
                                st0=st0, s=1,
                                precision = precision))
        simus$conf <- rnorm(n=nrow(simus),
                          mean=a * 0^(simus$response=="lower") +tau*simus$d,
                          sd=sqrt(tau))
      }
    } else {
      ### Simulation in the dynWEV model   ####
      w = paramDf$w
      svis = paramDf$svis
      if (all(SV==0)) {
        simus <- cbind(simus,
                     rdiffusion(n=nrow(simus), a=a, v=V[simus$condition]*simus$stimulus, t0 = t0,
                                z = z*a, d = 0, sz = a*sz, sv = 0,
                                st0=st0, s=1,
                                precision = precision))
        simus$evid_conf <- rnorm(n=nrow(simus),
                               mean=V[simus$condition]*simus$stimulus*tau*(-1)^(simus$response=="lower"),
                               sd = sqrt(tau))

      } else {
        simus$d <- rnorm(n=nrow(simus),
                       mean=V[simus$condition]*simus$stimulus,
                       sd = SV[simus$condition])
        simus <- cbind(simus,
                     rdiffusion(n=nrow(simus), a=a, v=simus$d, t0 = t0,
                                z = z*a, d = 0, sz = a*sz, sv = 0,
                                st0=st0, s=1,
                                precision = precision))
        simus$evid_conf <- rnorm(n=nrow(simus),
                               mean=simus$d*tau*(-1)^(simus$response=="lower"),
                               sd = sqrt(tau))
      }
      sigvis <- paramDf$sigvis
      simus$visibility <- rnorm(n=nrow(simus), mean= (simus$rt+tau)*abs(V[simus$condition]), sd = sqrt(svis^2*(tau+simus$rt)+(simus$rt+tau)^2*sigvis^2))
      simus$conf <- w*simus$evid_conf+(1-w)*simus$visibility
      simus$response <- if_else(simus$response=="upper", 1, -1)
    }
  } else {
    if (model =="2DSD") {
      w = -1
      svis = -1
      sigvis = -1
      vismu = rep(-1, nConds)
    } else {
      w = paramDf$w
      svis = paramDf$svis
      sigvis = paramDf$sigvis
      if ("vismu" %in% names(paramDf)) {
        vismu <- rep(paramDf$vismu, nConds)
      } else {
        vismu <- abs(V)
      }
    }

    simus <- expand.grid(condition = 1:nConds, stimulus=stimulus) %>%
      mutate(v  = V[.data$condition]*.data$stimulus,
             sv = SV[.data$condition],
             vismu = vismu[.data$condition]) %>%
      group_by(.data$condition, .data$stimulus) %>%
      summarise(as.data.frame(r_WEV(n=n, params=c(as.numeric(cur_data()[1:2]),a, z, sz, t0, st0, tau,w, as.numeric(cur_data()[3]), sigvis, svis),
                      model=which(model == c("2DSD", "dynWEV")),
                      delta = delta, maxT =maxrt), c("rt", "response", "conf"))) %>%
      rename(rt=3, response=4, conf=5)
  }

  ### Bin confidence measure for discrete ratings:
  symmetric_confidence_thresholds <- length(grep(pattern = "thetaUpper", names(paramDf), value = T))<1
  if (symmetric_confidence_thresholds) {
    nRatings <- length(grep(pattern = "^theta[0-9]", names(paramDf)))+1
  } else {
    nRatings <- length(grep(pattern = "^thetaUpper[0-9]", names(paramDf)))+1
  }
  if (model =="2DSD") {
    if (symmetric_confidence_thresholds) {
      thetas_upper <- c(-Inf, t(paramDf[,paste("theta",1:(nRatings-1), sep = "")]), Inf)
      thetas_lower <- a - thetas_upper
    } else {
      thetas_upper <- c(-Inf, t(paramDf[,paste("thetaUpper",1:(nRatings-1), sep = "")]), Inf)
      thetas_lower <- c(-Inf, t(paramDf[,paste("thetaLower",1:(nRatings-1), sep="")]), Inf)
    }
    levels_lower <- 6-cumsum(as.numeric(table(thetas_lower)))
    levels_lower <- levels_lower[-length(levels_lower)]
    levels_upper <- cumsum(as.numeric(table(thetas_upper)))
    levels_upper <- levels_upper[-length(levels_upper)]
    thetas_lower <- unique(thetas_lower)
    thetas_upper <- unique(thetas_upper)

    simus$rating <- 1
    simus$rating[simus$response==1] <- as.numeric(as.character(cut(simus$conf[simus$response==1],
                                                                          breaks=thetas_upper, labels = levels_upper)))
    simus$rating[simus$response==-1] <- as.numeric(as.character(cut(simus$conf[simus$response==-1],
                                                                             breaks=thetas_lower, labels=levels_lower)))
  } else {
    if (symmetric_confidence_thresholds) {
      thetas_upper <- c(-Inf, t(paramDf[,paste("theta",1:(nRatings-1), sep = "")]), Inf)
      thetas_lower <- c(-Inf, t(paramDf[,paste("theta",1:(nRatings-1), sep = "")]), Inf)
    } else {
      thetas_upper <- c(-Inf, t(paramDf[,paste("thetaUpper",1:(nRatings-1), sep = "")]), Inf)
      thetas_lower <- c(-Inf, t(paramDf[,paste("thetaLower",1:(nRatings-1), sep="")]), Inf)
    }
    levels_lower <- cumsum(as.numeric(table(thetas_lower)))
    levels_lower <- levels_lower[-length(levels_lower)]
    levels_upper <- cumsum(as.numeric(table(thetas_upper)))
    levels_upper <- levels_upper[-length(levels_upper)]
    thetas_lower <- unique(thetas_lower)
    thetas_upper <- unique(thetas_upper)

    simus$rating <- 1
    simus$rating[simus$response==1] <- as.numeric(as.character(cut(simus$conf[simus$response==1],
                                                              breaks=thetas_upper, labels = levels_upper)))
    simus$rating[simus$response==-1] <- as.numeric(as.character(cut(simus$conf[simus$response==-1],
                                                              breaks=thetas_lower, labels = levels_lower)))
  }
  simus$correct <- as.numeric(simus$response == simus$stimulus)
  if (simult_conf) {
    simus$rt <- simus$rt + tau
  }

  if (gamma==TRUE) {
    gamma_condition <- simus %>% group_by(.data$correct) %>%
      summarise(data.frame(t(rcorr.cens(.data$rating,S=.data$condition, outx=TRUE)))) %>%
      select(.data$correct, Gamma = .data$Dxy)
    gamma_rt <- simus %>% group_by(.data$correct) %>%
      summarise(data.frame(t(rcorr.cens(.data$rating,S=.data$rt, outx=TRUE))))%>%
      select(.data$correct, Gamma = .data$Dxy)
    gamma_correct <- simus %>% group_by(.data$condition) %>%
      summarise(data.frame(t(rcorr.cens(.data$rating,S=.data$correct, outx=TRUE))))%>%
      select(.data$condition, Gamma = .data$Dxy)
    gamma_rt_bycondition <- simus %>% group_by(.data$condition) %>%
      summarise(data.frame(t(rcorr.cens(.data$rating,S=.data$rt, outx=TRUE))))%>%
      select(.data$condition, Gamma = .data$Dxy)
    gamma_rt_byconditionbycorrect <- simus %>% group_by(.data$condition, .data$correct) %>%
      summarise(data.frame(t(rcorr.cens(.data$rating,S=.data$rt, outx=TRUE))))%>%
      select(.data$condition, Gamma = .data$Dxy)
  }
  if (agg_simus) {
    simus <- simus %>% group_by(.data$rating, .data$correct, .data$condition) %>%
      summarise(p = n()/(2*n)) %>%
      full_join(y=expand.grid(rating=1:nRatings, condition=1:nConds,
                              correct=c(0,1))) %>%
      mutate(p = ifelse(is.na(.data$p), 0, .data$p))
  } else {
    simus <- simus[c("condition", "stimulus", "response", "correct", "rt", "conf", "rating")]
  }
  if (gamma) {
    return(list("simus"=simus,
                "gamma" = list("condition" = gamma_condition,
                               "rt" = gamma_rt,
                               "correct" = gamma_correct,
                               "rt_bycondition" = gamma_rt_bycondition,
                               "rt_byconditionbycorrect" = gamma_rt_byconditionbycorrect)))
  } else {
    return(simus)
  }
}
