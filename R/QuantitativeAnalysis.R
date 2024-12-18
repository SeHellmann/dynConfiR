#' Quantitative Model Comparison
#'
#' \code{subject_modelweights} computes the model weights (as probabilities) of
#' individual subjects based on an information criterion (BIC, AIC, or AICc).
#' \code{group_BMS} performs a Bayesian model comparison based on marginal
#' likelihoods (alias model evidence), given for different models across different
#' subject on a group level using a fixed effects model and a random effects model
#' on the distribution of model probabilities (see Rigoux et al., 2014 and Details section).
#' \code{group_BMS_fits} is a wrapper for `group_BMS` that can be used with the
#' output of \code{\link{fitRTConfModels}}, i.e. a data frame with information
#' criteria for different models and subjects, using a information criterion to
#' approximate the model evidence.
#'
#' @param fits a data frame as returned by \code{\link{fitRTConfModels}}.
#' Should contain a column `model`indicating the model name, a column
#' `subject` (alternatively `sbj` or `participant`) indicating the grouping
#' structure of the data, and a column with the name given by the `measure`
#' argument containing the values of the information criterion that should be
#' used to approximate model evidence.
#' @param measure the name of the column indicating the information criterion
#' to approximate model evidence. For outputs of \code{\link{fitRTConfModels}},
#' the available measures are 'BIC', 'AIC', and 'AICc'. Any other approximation
#' for the model evidence may be used, the measure is transferred to log model
#' evidence by taking -measure/2.
#' @param opts a list with options for the iteration algorithm to estimate
#' the parameter of the Dirichlet distribution. Following values may be provided:
#' * \code{maxiter} the maximum number of iterations (Default: 200)
#' * \code{tol} the tolerance for changes in the free energy approximation
#'                 to stop the algorithm, if abs(FE(i+1)-FE(i))<tol the algorithm
#'                 is stopped (Default: 1e-4)
#' * \code{eps} The number to substitute values of 0 in calls to log (Default: 1e-32)
#' @param alpha0 a positive numeric vector representing the parameter of a
#' Dirichlet distribution used as prior over model probabilities. The length should
#' be equal to `nrow(mlp)` for `group_BMS`, and equal to the number of unique names
#' in the `model` column of `fits` for `group_BMS_fits`.
#' @param mlp a matrix containing the logarithm of marginal probabilities
#' (i.e. log model evidence) with N columns representing individuals (or any
#' other grouping structure) and K rows representing the models.
#'
#' @return \code{subject_modelweights} returns a data frame of subject-wise
#' model probabilities with rows for each subject and columns for the models
#' given by name and one column for the subject ID as given in the input.
#' \code{group_BMS} and \code{group_BMS_fits} return a list with two entries:
#'  * `model_weights`: a matrix with rows for each model (row names indicate the
#'                     model names for `group_BMS_fits` and for `group_BMS` if
#'                     row names are available in `mlp`), and following columns:
#'                     `alpha` (the alpha parameter of the Dirichlet posterior
#'                     over model probabilities in the population), `r` (the
#'                     mean probabilities of each model in the population), `ep`
#'                     and `pep` (exceedance and protected exceedance
#'                     probabilities for each model), and `fx_prob` (the
#'                     posterior model probabilities if a fixed true model is
#'                     assumed in the population).
#'  * `summary_stats`: a vector giving statistics for the Bayesian model comparison
#'                     that may be used for other analyses:
#'                     Bayesian omnibus risks: `bor` (random effects model against the
#'                     null model), `bor_fixed` (fixed effects model against the
#'                     null model), and `bor_re_fixed` (random effects model
#'                     against the fixed effects model), and
#'                     estimations of the Free Energy of the Dirichlet
#'                     distribution `FE` (random effects model), `FE0` (null model),
#'                    and `FEfixed` (fixed effects model)
#'
#' @details This set of function can be used for model comparisons on a group
#' level when the models were not fitted hierarchical but by fitting the models
#' independently to different subgroups (e.g. data from different subjects).
#'
#' The function `subject_modelweights` computes the model weights for each subject
#' separately to inspect predominant models but also heterogeneity within the
#' population.
#' The functions `group_BMS` and `group_BMS_fits` can be used for a Bayesian
#' model selection on the group level following the approach of Rigoux et al.
#' (2014).
#' The approach compares three different models for the generative structure of
#' the data and gives estimates for model probabilities for the fixed and random
#' effects models.

#' The **fixed effects model** assumes that there is a single model that
#' generated the data of all subjects. Thus, model weights may be computed directly
#' by multiplying the model weights computed on a subject-level. This model is
#' formulated in a Bayesian way using a Multinomial distribution over the models
#' as prior with some prior parameter alpha0 giving the prior model weights.
#' This is updated according to the marginal model likelihoods resulting in a
#' single poterior vector of model probabilities, reported in the column
#' `fx_prob` of the `model_weights` data frame.

#' The random effects model assumes that there is a vector of model probabilities
#' and each subject may generated by a different model, each drawn from a Multinomial
#' distribution. The Bayesian prior on this vector of model probabilities is
#' given by a Dirichlet distribution with some parameter alpha0. The function
#' uses a variational technique to approximate the alpha parameter of the
#' posterior Dirichlet distribution. Within this framework several statistics
#' may be used for model selection. The `model_weights` data frame reports the
#' posterior `alpha` parameter, as well as the posterior mean `r` of the corresponding
#' dirichlet distribution. The exceedance probability `ep` represents the probability
#' that given a random sample from the Dirichlet distribution the probability
#' of the model is greater than all other probailities. Finally, the
#' protected exceedance probability (`pep`) is a scaled version of the `ep`
#' multiplying the `ep` by one minus the Bayesian omnibus risk (BOR). The Bayesian omnibus
#' risk is the posterior probability of the  **null model** against the random
#' effects model. The **null model** assumes that all models are generating the
#' subjects' data with equal probability and results from taking the limit of
#' alpha0 towards infinity. The Bayesian omnibus risk is reported in the `summary_stats`
#' together with the free energy approximation of the null, the fixed effects,
#' and the random effects models.
#'
#' @references Rigoux, L., Stephan, K. E., Friston, K. J., & Daunizeau, J. (2014).
#'              Bayesian model selection for group studies - revisited. \emph{NeuroImage},
#'              84, 971–985. doi: 10.1016/j.neuroimage.2013.08.065
#'
# Daunizeau, J., Adam, V., & Rigoux, L. (2014). Vba: A probabilistic treatment
# of nonlinear models for neurobiological and behavioural data.
# \emph{PLOS Computational Biology}, 10(1), e1003441. doi: 10.1371/journal.pcbi.1003441
#'
#'
#' @author Sebastian Hellmann.
#'
#' @name QuantModelComparison
#' @aliases group_BMS, subject_modelweights, modelweights, group_BMS_fit, BMS
#'
#' @importFrom stats rgamma
#'
#' @examples
#' # Define a data frame with information criteria from model fits
#' # (this is a sub-data.frame from an output of fitRTConfModels with
#' # 8 subjects, three models and rounded information criteria)
#' fits <- data.frame(
#'   participant = rep(1:8, each=3),
#'   model = rep(c("dynaViTE", "2DSD", "PCRMt"), 8),
#'   BIC = c(5318, 5665, 1659, 3856, 5508, 3982, 3950, 3998,
#'           4114, 4216, 4314, 4419, 3170, 3489, 3256, 1950,
#'           1934, 2051, 3194, 3317, 3359, 9656, 10161, 4024),
#'   AIC = c(5211, 5577, 1577, 3750, 5420, 3899, 3843, 3911,
#'           4031, 4109, 4226, 4337, 3063, 3401, 3173, 1844,
#'           1847, 1969, 3087, 3229, 3277, 9549, 10074, 3942),
#'   AICc = c(5212, 5578, 1577, 3751, 5421, 3900, 3844, 3911,
#'            4032, 4110, 4227, 4337, 3064, 3402, 3174, 1845,
#'            1848, 1970, 3088, 3230, 3277, 9550, 10074, 3942))
#' # Compute subject-wise model probabitilities based on different ICs
#' subject_modelweights(fits, measure = "BIC")
#' subject_modelweights(fits, measure = "AIC")
#' subject_modelweights(fits, measure = "AICc")
#' # Conduct group-level Bayesian model selection based on BIC
#' group_BMS_fits(fits, measure="BIC")
#'
#'
#' ## General group-level Bayesian model selection based on any marginal log-probabilities
#' # Compute marginal log-likelihood based on BIC from fits
#' mlp <- matrix(NA, ncol=8, nrow=3)
#' for (i in 1:8) mlp[,i] <- fits[(i-1)*3 + 1:3, "BIC"]
#' mlp <- - mlp/(2)
#' rownames(mlp) <- c("dynaViTE", "2DSD", "PCRMt")
#' # conduct group BMS:
#' group_BMS(mlp)


### Compute subject-wise model weights using different information criteria
#' @rdname QuantModelComparison
#' @export
subject_modelweights <- function(fits, measure = "BIC"){
  if (!measure %in% c("BIC", "AIC", "AICc")) stop(paste0("measure=", measure, " must be present as column in fits"))
  models <- sort(unique(fits$model))
  sbj_col <- c("sbj", "subject", "participant")
  sbj_col <- sbj_col[which(sbj_col %in% names(fits))]
  participants <- unique(fits[[sbj_col]])
  n <- length(participants)
  K <- length(models)
  fits <- fits[order(fits[[sbj_col]], fits[["model"]]),]

  mlp <- matrix(NA, ncol=K, nrow=n)
  rownames(mlp) <- 1:n
  for (i in 1:n) {
    mlp[i,] <- fits[(i-1)*K + 1:K, measure]
    rownames(mlp)[i] <- unique(fits[(i-1)*K + 1:K, ][[sbj_col]])
  }
  colnames(mlp) <- fits$model[1:K]

  mlp <- - mlp/(2)
  ColMax <- apply(mlp, 1, max)
  mlp <- sweep(mlp, 1, ColMax)
  u_nk <- exp(mlp)
  #u_nk[is.infinite(u_nk)] <- 1e+64
  u_nk_sums <- rowSums(u_nk)

  model_weights <- as.data.frame(sweep(u_nk, 1, u_nk_sums, FUN="/"))
  model_weights[[sbj_col]] <- participants
  return(model_weights)
}


### Conduct Bayesian Group-Level model comparison
#' @rdname QuantModelComparison
#' @export
group_BMS_fits <- function(fits, measure = "BIC",
                             opts=list(), alpha0=NULL) {
  if (!measure %in% c("BIC", "AIC", "AICc")) stop(paste0("measure=", measure, " must be present as column in fits"))
  models <- sort(unique(fits$model))
  sbj_col <- c("sbj", "subject", "participant")
  sbj_col <- sbj_col[which(sbj_col %in% names(fits))]
  participants <- unique(fits[[sbj_col]])
  n <- length(participants)
  K <- length(models)
  fits <- fits[order(fits[[sbj_col]], fits[["model"]]),]
  mlp <- matrix(NA, ncol=n, nrow=K)
  for (i in 1:n) mlp[,i] <- fits[(i-1)*K + 1:K, measure]
  mlp <- - mlp/(2)
  rownames(mlp) <- models
  out <- group_BMS(mlp, opts, alpha0)
  return(out)
}

### Conduct Bayesian Group-Level model comparison
#' @rdname QuantModelComparison
#' @export
group_BMS <- function(mlp, opts=list(), alpha0=NULL) {
  if (is.null(opts$maxiter)) opts$maxiter <- 200 # max iterations for alpha-optimization
  if (is.null(opts$tol)) opts$tol <- 1e-4  # convergence criterion for alpha
  if (is.null(opts$eps)) opts$eps <- 1e-32   # replace 0s with eps in log-calls

  K <- nrow(mlp)
  n <- ncol(mlp)
  if (is.null(alpha0)) alpha0 <- rep(1, K)

  modelprobs_Null <- function(x) { # function for each n (x is a K-vector)
    g <- x - max(x)
    g <- exp(g)/sum(exp(g))
    res <- sum(g*(x - log(g+opts$eps) - log(K)))
    return(res)
  }
  FreeEnergyNull <- sum(apply(mlp, 2, modelprobs_Null))

  #% derive probabilities and free energy of the 'fixed-effect' analysis
  ss <-  rowSums(mlp) + log(1/K)
  logz <- ss - max(ss)
  z <-  exp(logz)/sum(exp(logz))
  fixed_effects_postprobs  <- z
  FreeEnergy_fixed <- sum(z * ss) - sum(z * log(z+opts$eps))

  alpha <- alpha0
  FreeEnergy <- 0
  #Full_FEs <- NULL
  #ln_p_y_mk # log(p(y_n | m_nk)) = - BIC/2
  for ( i in 1:opts$maxiter) {
    digamma_stuff <- digamma(alpha) - digamma(sum(alpha))
    logu_nk <- mlp + matrix(digamma_stuff, nrow=K, ncol=n, byrow=FALSE)
    ColMin <- apply(logu_nk, 2, min)
    logu_nk <- sweep(logu_nk, 2, ColMin)
    u_nk <- exp(logu_nk)

    u_nk[is.infinite(u_nk)] <- 1e+64
    model_sums <- colSums(u_nk)
    beta_nk <- sweep(u_nk, 2, model_sums, FUN="/")# z_nk bei Rigoux et al. (2014)
    # posterior.r in VBA toolbox

    Old_FreeEnergy <- FreeEnergy
    Sqf <- sum(lgamma(alpha)) - lgamma(sum(alpha)) - sum((alpha-1)*digamma_stuff)
    Sqm <- - sum(beta_nk * log(beta_nk + opts$eps))
    ELJ <- lgamma(sum(alpha0)) - sum(lgamma(alpha0)) + sum((alpha0-1)*digamma_stuff)
    ELJ <- ELJ + sum( beta_nk* (logu_nk))
    FreeEnergy <- Sqf + Sqm + ELJ

    beta_k <- rowSums(beta_nk)
    alpha <- alpha0 + beta_k
    #Full_FEs[i] <- FreeEnergy

    if (abs(Old_FreeEnergy-FreeEnergy)< opts$tol) { # equivalent to matlab (VBA-Toolbox TolFun = F(it)-F(it-1))
      # MaxIter in VBA toolbox is 32!
      #warning("Finished because of small update")
      break
    }
  }
  if (i == opts$maxiter) warning("Maximum number of iterations reached; alpha may have not converged!")

  ### Simulate exceedance probability
  simulate_ep <- function(n=100,alpha){
    # For dirichlet-rv we would scale Gamma-RV to a max of 1,
    # but we are interested in the index of the max only, so this it not
    # necessary
    K <- length(alpha)
    if (!is.null(names(alpha))) {
      res <- sapply(1:n, function(x) names(alpha)[which.max(rgamma(K, shape=alpha))])
      res <- factor(res, levels=names(alpha))
      res <- as.vector(table(res))/n
      names(res) <- names(alpha)
    } else {
      res <- sapply(1:n, function(x) which.max(rgamma(K, shape=alpha)))
      res <- as.vector(table(res))/n
    }
    return(res)
  }
  ep <- simulate_ep(1e+4, alpha)
  bor = 1/(1+exp(FreeEnergy-FreeEnergyNull));
  pep <- ep * (1 - bor) + bor / K

  bor_fixed <-     1/(1+exp(FreeEnergy_fixed-FreeEnergyNull));
  bor_re_fixed <-  1/(1+exp(FreeEnergy-FreeEnergy_fixed));

  model_weights <- matrix(c(alpha, alpha/sum(alpha), ep, pep, fixed_effects_postprobs), nrow=K)
  rownames(model_weights) <- rownames(mlp)
  colnames(model_weights) <- c("alpha", "r", "ep", "pep", "fx_prob")

  out <- list()
  out$model_weights <- model_weights
  out$summary_stats <- c(bor=bor, FE=FreeEnergy, FE0=FreeEnergyNull,
                         FEfixed=FreeEnergy_fixed, bor_fixed=bor_fixed, bor_re_fixed = bor_re_fixed)
  return(out)
}
