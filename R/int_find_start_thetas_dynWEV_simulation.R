get_thetas_for_init_grid_dynWEV_simulations <- function(init_grid, df, nRatings, simult_conf, fixed, mint0) {
  nConds <- length(unique(df$condition))
  conf_probs <- cumsum(table(df$rating))
  conf_probs <- conf_probs[1:(nRatings-1)]/conf_probs[nRatings]
  # df$correct <- as.numeric(df$stimulus == df$response)
  # MRT <- aggregate(rt~condition+correct, df, mean) %>%
  #   full_join(expand.grid(condition=levels(df$condition), correct=c(0,1))) %>%
  #   mutate(rt = ifelse(is.na(rt), 0, rt))
  # p_corrects <- aggregate(correct~condition, df, mean)[['correct']]

  if (!("t0" %in% names(fixed))) init_grid['t0'] <- init_grid['t0']*mint0


  get_start_thetas <- function(paramRow) {

    paramRow <- c(paramRow, unlist(fixed, use.names = TRUE))
    if (simult_conf & !("tau" %in% names(fixed))) {
      paramRow['tau'] <- paramRow['tau']*(mint0-paramRow['t0'])
    }
    V <- c(t(paramRow[paste("v", 1:nConds, sep="")]))
    conf <- with(as.data.frame(as.list(paramRow)), rdynaViTE(800, a, c(V,-V), t0, z, 0, sz, sv, st0, tau, w, abs(c(V,V)), sigvis, svis, lambda, 1))$conf
    thetas <- quantile(conf, probs=conf_probs, names = FALSE)
    c(thetas[1],diff(thetas))
  }


  init_thetas <- apply(init_grid, FUN=get_start_thetas, MARGIN=1) # , simplify = TRUE
  init_thetas <- t(init_thetas)
  init_thetas[,2:(nRatings-1)] <- pmax(init_thetas[,2:(nRatings-1)], 0.001)
  init_thetas

}
