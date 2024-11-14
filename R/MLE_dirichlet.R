#' Estimate the parameter of a Dirichlet distribution
#'
#' The function `MLE_dirichlet` performs a maximum-likelihood estimation of the
#' \eqn{\alpha}{\alpha} parameter of a Dirichlet distribution for a given sample of
#' probability vectors.
#'
#' @param probs a matrix with N rows representing observations of probability
#' vectors and K columns representing the classes. Therefore, values of each row
#' should sum to 1.
#' @param alpha0 vector of K=ncol(probs) values as starting parameter for the optimization.
#' Values have to be greater 0.
#'
#' @return Returns a numeric vector of length K=ncol(probs) representing the
#' \eqn{\alpha}{\alpha} of the Dirichlet distribution.
#'
#' @details The density of the Dirichlet distribution for
#' \eqn{\alpha = (\alpha_1, ..., \alpha_K )} and
#' \eqn{\alpha_i > 0 \forall i=1,...,K} is given by
#' \deqn{f(p|\alpha)=\frac{1}{B(\alpha)} \prod_{i=1}{K} p_{i}^{\alpha_i - 1},}
#' if \eqn{0\leq p_i \leq 1 \forall i = 1,...,K} and \eqn{\sum_{i=1}^{K} p_i ) 1},
#' and \eqn{f(p|\alpha) = 0}, else.
#'
#' The function optimizes the log-likelihood of a sample of probability vectors
#' given in `probs` using the function \code{\link{optim}} and a Nelder-Mead
#' algorithm.
#'
#' @author Sebastian Hellmann.
#'
#' @name MLE_dirichlet
#' @importFrom stats optim
#'
#' @examples
#' probs <- matrix(c(0.2, 0.4, 0.2, 0.4, 0, 0.2, 0.4, 0.4, 0.6, 0.2, 0.2,
#'                   0.4, 0.4, 0.2, 0.2, 0.4, 0.8, 0.4), ncol=3)
#' MLE_dirichlet(probs)


### Compute participant-wise model weights using different information criteria
#' @rdname MLE_dirichlet
#' @export
MLE_dirichlet <- function(probs, alpha0=rep(1, ncol(probs))) {
  #probs <- as.matrix(model_weights[,-ncol(model_weights)])
  K <- ncol(probs)
  N <- nrow(probs)
  probs[probs ==0] <- .Machine$double.xmin
  probs <- log(probs)
  #probs[probs == -Inf] <- -.Machine$double.xmax
  log_p_k <- apply(probs, 2, "mean")
  log_likelikood_unbound <- function(logalpha) {
    # N*
    -1*(lgamma(sum(exp(logalpha)))- sum(lgamma(exp(logalpha))) + sum((exp(logalpha)-1)*log_p_k))
  }
  optim_free <- optim(log(alpha0), log_likelikood_unbound)
  return(exp(optim_free$par))
}
