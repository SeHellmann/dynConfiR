get_thetas_for_init_grid_dynWEV <- function(init_grid, df, nRatings, simult_conf, fixed) {
  nConds <- length(unique(df$condition))
  conf_probs <- cumsum(table(df$rating))
  conf_probs <- conf_probs[1:(nRatings-1)]/conf_probs[nRatings]
  df$correct <- as.numeric(df$stimulus == df$response)
  MRT <- aggregate(rt~condition+correct, df, mean) %>%
    full_join(expand.grid(condition=levels(df$condition), correct=c(0,1))) %>%
    mutate(rt = ifelse(is.na(rt), 0, rt))
  p_corrects <- aggregate(correct~condition, df, mean)[['correct']]

  if (!("t0" %in% names(fixed))) init_grid['t0'] <- init_grid['t0']*mint0


  get_start_thetas <- function(paramRow) {

    paramRow <- c(paramRow, unlist(fixed, use.names = TRUE))
    if (simult_conf & !("tau" %in% names(fixed))) {
      paramRow['tau'] <- paramRow['tau']*(mint0-paramRow['t0'])
    }
    V <- c(t(paramRow[paste("v", 1:nConds, sep="")]))
    MRT_corr <- pmax(0,(filter(MRT, .data$correct==1)$rt-paramRow['t0']-paramRow['st0']/2))
    MRT_false <- pmax(0,(filter(MRT, .data$correct==0)$rt-paramRow['t0']-paramRow['st0']/2))

    if (simult_conf) {
      MRT_corr <-  pmax(0,MRT_corr -paramRow['tau'])
      MRT_false <- pmax(0,MRT_false -paramRow['tau'])
    }

    MRT_corr_ttau <- MRT_corr + paramRow[['tau']]
    MRT_false_ttau <- MRT_false + paramRow[['tau']]


    Mconf_false = (-paramRow[['w']]*
                     (paramRow[['tau']]*V-paramRow[['a']]*paramRow[['z']]*(paramRow[['sv']]^2*MRT_false_ttau+1))/(1+paramRow[['sv']]^2*MRT_false)+
                     (1-paramRow[['w']])*MRT_false_ttau*V) / MRT_false_ttau^paramRow[['omega']]

    Mconf_corr = (-paramRow[['w']]*
                     (-paramRow[['tau']]*V-paramRow[['a']]*(1-paramRow[['z']])*(paramRow[['sv']]^2*MRT_corr_ttau+1))/(1+paramRow[['sv']]^2*MRT_corr)+
                     (1-paramRow[['w']])*MRT_corr_ttau*V) / MRT_corr_ttau^paramRow[['omega']]

    VRconf_false = (paramRow[['w']]^2*
                      (paramRow[['tau']] +
                         (paramRow[['sv']]^2*paramRow[['tau']]^2)/
                         (1+paramRow[['sv']]^2*MRT_false)) +
                      ((1-paramRow[['w']])^2*
                         (paramRow[['svis']]^2*MRT_false_ttau + MRT_false_ttau^2 * paramRow[['sigvis']]^2)))/MRT_false_ttau^(2*paramRow[['omega']])

    VRconf_corr = (paramRow[['w']]^2*
                      (paramRow[['tau']] +
                         (paramRow[['sv']]^2*paramRow[['tau']]^2)/
                         (1+paramRow[['sv']]^2*MRT_corr)) +
                      ((1-paramRow[['w']])^2*
                         (paramRow[['svis']]^2*MRT_corr_ttau + MRT_corr_ttau^2 * paramRow[['sigvis']]^2)))/MRT_corr_ttau^(2*paramRow[['omega']])
    if (any(is.na(Mconf_corr))||any(is.na(Mconf_false))||any(is.na(VRconf_corr))||any(is.na(VRconf_false))) {
      print(paramRow)
      stop("Some values are NA")
    }
    mixcdf <- function(conf) 1/nConds * sum((1-p_corrects)* pnorm(conf, mean=Mconf_false, sd=sqrt(VRconf_false))+
                                                       p_corrects*pnorm(conf, mean=Mconf_corr, sd=sqrt(VRconf_corr)))
    thetas <- NULL
    for (i in 1:length(conf_probs)) {
      thetas[i] <- optimize(function(conf) (mixcdf(conf)-conf_probs[i])^2,
                            lower=min(Mconf_false,Mconf_corr)- 4*max(c(VRconf_false, VRconf_corr)),
                            upper=max(Mconf_false,Mconf_corr)+ 4*max(c(VRconf_false, VRconf_corr)))$minimum
    }
    c(thetas[1],diff(thetas))
  }


  init_thetas <- apply(init_grid, FUN=get_start_thetas, MARGIN=1) # , simplify = TRUE
  init_thetas <- t(init_thetas)
  init_thetas[,2:(nRatings-1)] <- pmax(init_thetas[,2:(nRatings-1)], 0.001)
  init_thetas

}
