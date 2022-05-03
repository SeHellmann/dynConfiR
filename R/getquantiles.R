#' Get Quantiles from predicted Density or CDF
#'
#' \code{compute_quantiles_cdf} computes quantiles for a given CDF vector.
#' \code{RTDensityToQuantiles} computes the quantiles from a vector of probability
#' density values within groups of other given variables (in the context of confidence
#' models: conditions, correct/incorrect answers and confidence ratings), if
#' available.
#'
#' @param pred dataframe. Should have at least two columns: rt (for reaction times)
#' and dens or densscaled. All other columns will be used as grouping factors.
#' As input for pred, the outputs of predictRT and predictRTModels can be used.
#' @param p numeric vector. Probabilities for returned quantiles. Default:
#' c(.1, .3, .5, .7, .9).
#' @param agg_over character. Names of columns to aggregate over (using the mean of
#' densities, which is valid only, if groups occur with equal probabilities) before
#' computing the quantiles.
#' @param cdf numeric. A increasing vector of the same length as \code{rt} giving the CDF for RT-Values.
#' @param rt numeric. A increasing vector of same length as \code{cdf}. May be anything, not only rt's.

#'
#' @return \code{RTDensityToQuantiles} returns a tibble with columns p and q indicating
#' probabilities and respective quantiles. Furhtermore, the output has grouping columns
#' identical to the additional columns in the input (without rt, dens and densscaled),
#' but without the ones in the agg_over argument. \code{compute_quantiles_cdf}
#' returns only a data.frame with columns p and q.
#'
#' @details The RTs in the input data frame should be equidistant, i.e. the difference
#' between consecutive RT-values should be constant. For a reasonable accuracy the number
#' of steps in the input should be very high, i.e. the distant between rt values small.
#' The precision of the output is bounded from above by the step size of the RTs in the input.
#'
#' If available, the column densscaled will be preferred.
#' If the column dens is used, it will be scaled by the summed probability.
#'
#' Attention should be given to the columns of pred other then rt and dens/densscaled, because
#' the quantiles are computed
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
#' @name RTDensityToQuantiles
#' @import dplyr
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @aliases getquantiles getRTquantiles
#' @importFrom Rcpp evalCpp
#'


#' @rdname RTDensityToQuantiles
#' @export
RTDensityToQuantiles <- function(pred, p = c(.1,.3,.5,.7,.9),
                                 agg_over = NULL){
  pred <- ungroup(pred)



  if ("densscaled" %in% names(pred)) {
    if ("dens"  %in% names(pred)) {
      pred <- select(pred, -"dens")
    }
    temp <- pred %>% select(-rt, -densscaled) %>%
      group_by(across()) %>% summarise(N=n())
    if (min(temp$N) < 100) {
      warning(paste("There are only", min(temp$N), "rows for at least one subgroup of the data set.",
      "\nConsider refining the rt-grid for more accurate computations."))
    }
    pred <- pred %>% group_by(pred[,setdiff(names(pred), c("rt", "densscaled"))]) %>%
      arrange(.data$rt) %>%
      mutate(dt = (c(0, diff(.data$rt))+c(diff(.data$rt),0))*0.5)
  } else {
    if (!("dens" %in% names(pred))) {
      stop("At least dens or densscaled should be a column of pred")
    }
    temp <- pred %>% select(-rt, -dens) %>%
      group_by(across()) %>% summarise(N=n())
    if (min(temp$N) < 100) {
      warning(paste("There are only", min(temp$N), "rows for at least one subgroup of the data set.",
                    "\nConsider refining the rt-grid for more accurate computations."))
    }

    pred <- pred %>% group_by(pred[,setdiff(names(pred), c("rt", "dens"))]) %>%
      arrange(.data$rt) %>%
      mutate(dt = (c(0, diff(.data$rt))+c(diff(.data$rt),0))*0.5,
             p=sum(.data$dens*.data$dt),
             densscaled = .data$dens/.data$p) %>%
      select(-c("p", "dens"))
  }

  pred <- pred %>% group_by(pred[, setdiff(names(pred), c("rt", "densscaled", "dt"))]) %>%
    mutate(cdfscaled= cumsum(.data$densscaled*.data$dt)) %>%
    select(-"densscaled", - "dt")
  if (!is.null(agg_over)) {
    pred <- pred %>% group_by(pred[,setdiff(names(pred), c("cdfscaled", agg_over))]) %>%
      summarise(cdfscaled = mean(.data$cdfscaled))
  }
  pred <- pred %>% group_by(pred[, setdiff(names(pred), c("rt", "cdfscaled"))]) %>%
    summarise(compute_quantiles_cdf(.data$cdfscaled, .data$rt, p = p))

  pred
}


#' @rdname RTDensityToQuantiles
#' @export
compute_quantiles_cdf <- function(cdf, rt, p) {
  q <- rep(0, length(p))
  for (i in 1:length(p)) {
    loc <- min(which(cdf>= p[i]))
    q[i] <- rt[loc]
  }
  return(data.frame(p = p, q = q))
}
