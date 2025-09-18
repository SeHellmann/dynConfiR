get_thetas_for_init_grid_MTLNR_simulations <- function(init_grid, df, nRatings, fixed) {
  nConds <- length(unique(df$condition))
  conf_probs <- cumsum(table(df$rating))
  conf_probs <- conf_probs[1:(nRatings-1)]/conf_probs[nRatings]


  # paramRow <- as.numeric(init_grid[1,])
  # names(paramRow) <- names(init_grid)
  get_start_thetas <- function(paramRow) {

    paramRow <- c(paramRow, unlist(fixed, use.names = TRUE))
    V <- c(t(paramRow[paste("v", 1:nConds, sep="")]))
    # if (paramRow["a"] == "b") paramRow["a"] <- paramRow["b"]
    # if (paramRow["b"] == "a") paramRow["b"] <- paramRow["a"]
    conf <- with(as.data.frame(as.list(paramRow)), rMTLNR(800, c(5,5), # Arbitrary thresholds
                                                          c(V,rep(0, nConds)), # mu_v1 for stim 1 and 2
                                                          c(rep(0, nConds),V), # mu_v2 for stim 1 and 2
                                                          s_v1, s_v2, rho_v, mu_d1, mu_d2, s_d1, s_d2, rho_d,
                                                          t0=0, st0=0))$conf
    thetas <- quantile(conf, probs=conf_probs, names = FALSE)
    c(thetas[1],diff(thetas))
  }


  init_thetas <- apply(init_grid, FUN=get_start_thetas, MARGIN=1) # , simplify = TRUE
  init_thetas <- t(init_thetas)
  init_thetas[,1] <- pmax(init_thetas[,1], 1.001)
  init_thetas[,2:(nRatings-1)] <- pmax(init_thetas[,2:(nRatings-1)], 0.001)
  init_thetas

}
