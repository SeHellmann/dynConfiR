compute_parammatrix <- function(beta, model_matrix, collapse=FALSE) {
  #parnames <- c("v", "z", "a", "d", "sz", "t0", "st0", "sv", "tau", "w", "svis", "sigvis", "lambda", "s")
  parnames <- c("v", "z", "a", "sz", "t0", "st0", "d", "sv", "tau", "w", "svis", "sigvis", "lambda", "s", "th1", "th2")
  parammatrix <-  matrix(NA, nrow=nrow(model_matrix), ncol=length(parnames))
  colnames(parammatrix) <- parnames
  fixed01 <- c("z", "sz", "w", "d")
  fixedpos <- c("v", "a", "t0", "st0", "sv", "svis", "sigvis", "lambda", "s")

  beta_names <- names(beta)

  for (i in 1:length(parnames)) {
    cur_par <- parnames[i]
    beta_par <- grep(pattern=paste0(cur_par, "_"), beta_names, value=TRUE)
    if (length(beta_par)>0) {

      parammatrix[,cur_par] <-
        model_matrix[,sub(paste0(cur_par, "_"), beta_par), drop=FALSE] %*%
        beta[beta_par]
    } else {
      parammatrix[,cur_par] <- beta[cur_par]
    }
    if (cur_par %in% fixed01) parammatrix[,cur_par] <- pnorm(parammatrix[,cur_par])
    if (cur_par %in% fixedpos) parammatrix[,cur_par] <- exp(parammatrix[,cur_par])
    # if (cur_par == "sz") parammatrix[,"sz"] <- (pmin(parammatrix[,"z"], (1-parammatrix[,"z"]))*2)*parammatrix[,"sz"]
    # if (cur_par == "t0") parammatrix[,"t0"] <- maxt0*parammatrix[,"t0"]
    # if (cur_par == "d") parammatrix[,"d"] <- parammatrix[,"t0"]*parammatrix[,"d"]
    # if (cur_par == "tau") {
    #   if (restr_tau == Inf) {
    #     parammatrix[,cur_par] <- exp(parammatrix[,cur_par])
    #   } else if (simult_conf) {
    #     parammatrix[,cur_par] <- pnorm(parammatrix[,cur_par])*(maxt0-parammatrix[,"t0"])
    #   } else {
    #     parammatrix[,cur_par] <- restr_tau * pnorm(parammatrix[,cur_par])
    #   }
    # }
  }
  #parammatrix[,"v"] <- parammatrix[,"v"]*DVs$stimulus

  par_thetas <- grep(names(beta), pattern="theta", value=TRUE)
  if (sym_thetas) {
    Thetas <- beta[par_thetas]
    Thetas <- c(-1e+21, Thetas[1] + c(0, cumsum(exp(Thetas[-1]))), Inf)
    Thetas <- c(Thetas, Thetas)
  } else {
    ThetasLower <- beta[grep(par_thetas, pattern="Lower", value = TRUE)]
    ThetasLower <- c(-1e+21, ThetasLower[1] + c(0, cumsum(exp(ThetasLower[-1]))), 1e+21)
    ThetasUpper <- beta[grep(par_thetas, pattern="Upper", value = TRUE)]
    ThetasUpper <- c(-1e+21, ThetasUpper[1] + c(0, cumsum(exp(ThetasUpper[-1]))), 1e+21)
    Thetas <- c(ThetasLower, ThetasUpper)
  }
  parammatrix[,c("th1", "th2")] <- cbind(Thetas[(DVs$response==1)*(nRatings+1) + DVs$rating],
                                         Thetas[(DVs$response==1)*(nRatings+1) + DVs$rating+1])
}
