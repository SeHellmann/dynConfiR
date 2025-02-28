get_thetas_for_init_grid_RMs_simulations <- function(init_grid, df, model, nRatings, fixed, time_scaled, fitted_weights) {
  nConds <- length(unique(df$condition))
  conf_probs <- cumsum(table(df$rating))
  conf_probs <- conf_probs[1:(nRatings-1)]/conf_probs[nRatings]


  # paramRow <- as.numeric(init_grid[1,])
  # names(paramRow) <- names(init_grid)
  get_start_thetas <- function(paramRow) {

    paramRow <- c(paramRow, unlist(fixed, use.names = TRUE))
    V <- c(t(paramRow[paste("v", 1:nConds, sep="")]))
    if (paramRow["a"] == "b") paramRow["a"] <- paramRow["b"]
    if (paramRow["b"] == "a") paramRow["b"] <- paramRow["a"]
    if (!time_scaled) {
      paramRow["wx"] <- 1
      paramRow["wrt"] <- 0
      paramRow["wint"] <- 0
    } else {
      if (length(fitted_weights) == 1) {
        paramRow[fitted_weights] <- paramRow[fitted_weights]*(1-fixed[[grep("^w", names(fixed), value = TRUE)]])
        paramRow[grep("^w", names(fixed), value = TRUE)] <- fixed[[grep("^w", names(fixed), value = TRUE)]]
        paramRow[setdiff(c("wx", "wrt", "wint"),names(paramRow))] <- 1- sum(as.numeric(paramRow[grep("^w", names(paramRow), value=TRUE)]))
      }
      if (length(fitted_weights)==2) {
        paramRow["wrt"] <- paramRow["wrt"]*(1-paramRow["wx"])
        paramRow["wint"] <- 1 - paramRow["wx"] - paramRow["wrt"]
      }
    }

    if (model=="IRM") {
      conf <- with(as.data.frame(as.list(paramRow)), rIRM(800, c(V,-V),c(-V,V), a, b, wx, wrt, wint, t0=0, st0=0,
                                                          time_scaled=time_scaled))$conf
    } else {
      conf <- with(as.data.frame(as.list(paramRow)), rPCRM(800, c(V,-V),c(-V,V), a, b, wx, wrt, wint, t0=0, st0=0,
                                                           time_scaled=time_scaled))$conf
    }
#    conf <- with(as.data.frame(as.list(paramRow)), rdynaViTE(800, a, c(V,-V), t0, z, 0, sz, sv, st0, tau, w, abs(c(V,V)), sigvis, svis, lambda, 1))$conf
    thetas <- quantile(conf, probs=conf_probs, names = FALSE)
    c(thetas[1],diff(thetas))
  }


  init_thetas <- apply(init_grid, FUN=get_start_thetas, MARGIN=1) # , simplify = TRUE
  init_thetas <- t(init_thetas)
  init_thetas[,2:(nRatings-1)] <- pmax(init_thetas[,2:(nRatings-1)], 0.001)
  init_thetas

}
