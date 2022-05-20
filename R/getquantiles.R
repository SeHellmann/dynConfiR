#' Get Quantiles from vectors of PDF or CDF values
#'
#' `CDFtoQuantiles` computes quantiles for a given CDF.
#' `PDFtoQuantiles` computes the quantiles for given PDF values within
#' groups of other variables, if available.

#' @param pdf_df dataframe. Should have at least two columns:
#' * `rt` (for reaction times) or `x` for the support values of the pdf
#' * `dens` or `pdf` for the pdf values
#' * All other columns will be used as grouping factors, for which separate quantiles will be returned.
#' @param p numeric vector. Probabilities for returned quantiles. Default:
#' c(.1, .3, .5, .7, .9).
#' @param agg_over character. Names of columns in `pdf_df` to aggregate over (using the mean of
#' densities, which is valid only, if groups occur with equal probabilities) before
#' computing the quantiles.
#' @param scaled logical. Indicating whether the pdf values are from a proper probability
#' distribution. Non-scaled pdfs will scaled to 1. If `scaled` is TRUE, this may cause
#' problems with high probabilities. In any case we strongly recommend to cover the most
#' probability mass with the values in the support vector.
#' @param cdf numeric. A increasing vector of the same length as `x` giving the CDF for respective x-Values.
#' Dataframe inputs are accepted. If a column `x` is available there, this will be used as support values.
#' @param x numeric. A increasing vector of same length as `cdf`. Can also be specified as column of `cdf`.
#'
#' @return `PDFtoQuantiles` returns a tibble with columns p and q indicating
#' probabilities and respective quantiles. Furthermore, the output has grouping columns
#' identical to the additional columns in the input (without rt/x, dens and densscaled),
#' but without the ones in the agg_over argument. `CDFtoQuantiles`
#' returns only a data.frame with columns p and q.
#'
#' @details
#' For a reasonable accuracy the number of steps in the support column (`rt`/`x`)
#' should be high, i.e. the distance between values small.
#' We recommend, to ensure that the support vector in the input to be equidistant,
#' i.e. the difference between consecutive support values should be constant, though
#' this is not required.
#' If both column names `x` and `rt` are present in `pdf_df`, `rt` will be preferred.
#' Attention should be given to the columns of `pdf_df` other than `rt`/`x`
#' and `dens`/`pdf`.
#'
#' The column for the pdf may be scaled to integrate to 1 but do not have to.
#'
#' ## Quantile computation in the `dynConfiR` package
#' As argument `pdf_df`, the outputs of `predictRT` and `predictRTModels` from the
#' `dynConfiR` package can be used. In the context of confidence models grouping factors
#' often used are conditions, correct/incorrect answers and confidence ratings.
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
#' @name PDFtoQuantiles
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @import dplyr
#'
#' @aliases getquantiles getRTquantiles RTDensityToQuantiles
#' @importFrom Rcpp evalCpp
#'
#' @examples
#' ## Demonstrate PDFtoQuantiles
#' pred <- expand.grid(model = c("dynWEV", "PCRMt"),
#'                     rt =  seq(0, 15, length.out=1200),
#'                     condition = c(1,2,3),
#'                     rating = c(1,2))
#' pred$dens <- dchisq(pred$rt, 3) # pdf may also be used as column name
#' head(pred)
#' res <- PDFtoQuantiles(pred, p=c(0.3, 0.5, 0.7))
#' head(res)
#' nrow(res) #= 3(quantiles)*2(models)*3(conditions)*2(rating)
#' # Compare to true quantiles of Chi-square distribution
#' qchisq(p=c(0.3, 0.5, 0.7), 3)
#' res$q[1:3]
#'
#'
#' res2 <- PDFtoQuantiles(pred, p=c(0.3, 0.5, 0.7), agg_over = "model")
#' nrow(res2) #=18 because res aggregated over models
#'
#' \dontrun{
#'   pred$pdf <- dchisq(pred$rt, 3)
#'   head(pred)
#'   # following call throws a warning, because both columns pdf and dens are present
#'   PDFtoQuantiles(pred, p=c(0.3, 0.5, 0.7), agg_over = "model")
#' }
#'
#' \dontrun{
#'   pred2 <- data.frame(rt=seq(0, 7, length.out=100))
#'   pred2$dens <- dchisq(pred2$rt, 5)
#'   # following gives a warning, because density is assumed to be scaled (scaled=TRUE), i.e.
#'   # integrate to 1, but the .95 quantile is not reached in the rt column
#'   PDFtoQuantiles(pred2, p=c(0.3, 0.5, 0.95), scaled=TRUE) # Gives a warning
#' }
#'
#' ## Demonstrate CDFtoQuantiles
#' X <- seq(-2, 2, length.out=300)
#' pdf_values <- pnorm(X)
#' CDFtoQuantiles(pdf_values, X, p=c(0.2, 0.5, 0.8))
#' qnorm(c(0.2, 0.5, 0.8))





#' @rdname PDFtoQuantiles
#' @export
PDFtoQuantiles <- function(pdf_df, p = c(.1,.3,.5,.7,.9),
                                 agg_over = NULL, scaled=FALSE){
  pdf_df <- ungroup(pdf_df)
  columns <- names(pdf_df)
  if ("x" %in% columns) {
    if ("rt" %in% columns) {
      warning("The column 'rt' is used as support values, not 'x'!")
    } else {
      pdf_df$rt <- pdf_df$x
    }
  } else if (!("rt" %in% columns)) stop("Either 'x' or 'rt' must be a column of pdf_df")

  if ("pdf" %in% columns) {
    if ("dens" %in% columns) {
      warning("The column 'dens' is used as density values, not 'pdf'!")
      pdf_df$pdf <- NULL
    } else {
      pdf_df$dens <- pdf_df$pdf
      pdf_df$pdf <- NULL
    }
  } else if (!("dens" %in% columns)) stop("Either 'dens' or 'pdf' must be a column of pdf_df")
  if ("densscaled" %in% names(pdf_df)) pdf_df$densscaled <- NULL

  temp <- pdf_df %>% select(-c("rt", "dens")) %>%
    group_by(across()) %>% summarise(N=n())
  if (min(temp$N) < 100) {
    warning(paste("There are only", min(temp$N), "rows for at least one subgroup of the data set.",
                  "\nConsider refining the rt-grid for more accurate computations."))
  }
  if (!scaled) {
    pdf_df <- pdf_df %>% group_by(pdf_df[,setdiff(names(pdf_df), c("rt", "dens"))]) %>%
      arrange(.data$rt) %>%
      mutate(dt = (c(0, diff(.data$rt))+c(diff(.data$rt),0))*0.5,
             p=sum(.data$dens*.data$dt),
             densscaled = .data$dens/.data$p) %>%
      select(-c("p", "dens"))
  } else {
    pdf_df <- pdf_df %>% group_by(pdf_df[,setdiff(names(pdf_df), c("rt", "dens"))]) %>%
      arrange(.data$rt) %>%
      mutate(dt = (c(0, diff(.data$rt))+c(diff(.data$rt),0))*0.5,
             densscaled = .data$dens) %>%
      select(-c("dens"))
  }


  pdf_df <- pdf_df %>% group_by(pdf_df[, setdiff(names(pdf_df), c("rt", "densscaled", "dt"))]) %>%
    mutate(cdfscaled= cumsum(.data$densscaled*.data$dt)) %>%
    select(-"densscaled", - "dt")
  if (!is.null(agg_over)) {
    pdf_df <- pdf_df %>% group_by(pdf_df[,setdiff(names(pdf_df), c("cdfscaled", agg_over))]) %>%
      summarise(cdfscaled = mean(.data$cdfscaled))
  }
  pdf_df <- pdf_df %>% group_by(pdf_df[, setdiff(names(pdf_df), c("rt", "cdfscaled"))]) %>%
    summarise(CDFtoQuantiles(.data$cdfscaled, .data$rt, p = p))

  pdf_df
}


#' @rdname PDFtoQuantiles
#' @export
CDFtoQuantiles <- function(cdf, x=NULL, p) {
  if (is.data.frame(cdf)) {
    if (!("cdf" %in% names(cdf))) stop("cdf is a data.frame but neither 'cdf' nor 'p' is a column")
    if (ncol(cdf)==2) {
      if ("p" %in% names(cdf)) {
        x <- cdf[,-which(names(cdf)=="p")]
        cdf <- cdf[["p"]]
      } else {
        x <- cdf[,-which(names(cdf)=="cdf")]
        cdf <- cdf[["cdf"]]
      }
    } else if (ncol(cdf)>2) {
      if (!("x" %in% names(cdf))) stop("cdf is a data.frame with more than two columns but 'x' is not a column. \n Try using the function with two vector arguments.")
      x <- cdf[['x']]
      if ("p" %in% names(cdf)) {
        cdf <- cdf[["p"]]
      } else {
        cdf <- cdf[["cdf"]]
      }
    }
  } else {
    if (is.null(x)) stop("'x' must be specified unless 'cdf' is a data.frame")
  }
  q <- rep(0, length(p))
  for (i in 1:length(p)) {
    loc <- min(which(cdf>= p[i]))
    q[i] <- x[loc]
  }
  return(data.frame(p = p, q = q))
}
