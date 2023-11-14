## Confidence thresholds (theta) are fitted as decision time thresholds!!!
## (i.e. increasing "confidence" here, is increasing RT and thus corresponds to
## lower confidence)

fittingDDMConf<- function(df, nConds, nRatings, fixed, sym_thetas,
                        grid_search, init_grid=NULL, opts,
                        logging, filename,
                        useparallel, n.cores, precision,
                        used_cats, actual_nRatings,st0stepsize){
  ## Be sure that the parallel cluster is stopped if anything happens (error or user interrupt)
  on.exit(try(stopCluster(cl), silent = TRUE))
  ## Set restrictions on st0
  mint0 <-  min(df$rt)
  #if (mint0 < 0.2) {mint0 <- 0.2}
  common_RTrange <- df %>% group_by(.data$rating) %>%
    summarise(minRT = min(.data$rt), maxRT=max(.data$rt))
  common_RTrange <- c(sort(common_RTrange$maxRT, decreasing = TRUE)[2],
    sort(common_RTrange$minRT, decreasing = FALSE)[2])
  minst0 <- common_RTrange[1]- common_RTrange[2]
  st0 <- c(minst0, (3*minst0+max(df$rt))/4, (minst0+max(df$rt))/2, (minst0+4*max(df$rt))/5)
  ## Guess suitable confidence thresholds from theoretical distribution of
  ## the confidence measure and proportion of ratings in the data

  # observed rt are equal to dt+nondectime, and all dts from rating=1 have to be
  # dt > thetaMax --> rt = dt + nondectime  > thetaMax + nondectime > thetaMax
  # --> therefore: thetaMax <= min(rt | rating=1) - t0
  thetaMax <- min(filter(df, .data$rating == 1)$rt)
  thetaMax <- seq(0.6, 0.9, length.out=3)*thetaMax
  # observed rt are maximal equal to dt+t0+st0, and all dts from rating=5 have to be
  # faster than thetaMin (i.e. dt < thetaMin)
  # --> rt = dt + nondectime < thetaMin + nondectime < thetaMin + t0 + st0
  # --> therefore: thetaMin > rt - t0 - st0 for all rts from rating=5
  # --> therefore: thetaMin > max(rt | rating=5) - t0 - st0
  if (nRatings > 2) {
    # theta_(nRatings-1) = thetaMin * thetaMax
    thetaMin <- c(1.02, 1.2, 1.4) * max(filter(df, .data$rating==5)$rt)
  } else {
    thetaMin <- 1.1 * max(filter(df, .data$rating == 5)$rt)
  }

  ### 1. Generate initial grid for grid search over possible parameter sets ####
  #### Create grid ####
  if (is.null(init_grid)) {
    init_grid <- expand.grid(a = c(0.8,  1.2, 2, 3),
                             vmin = c(0.01, 0.5, 1.3),
                             vmax = c(2, 4, 6),
                             sv = c(0.1, 1.5),
                             z = c(0.4, 0.6),
                             sz = c(0.01, 0.1, 0.3),
                             t0 = c(0.01, max(mint0-1.3, 0.02), max(min(mint0-1,0.2),mint0/2)),
                             st0 = st0,
                             thetaMax = thetaMax,
                             thetaMin = thetaMin)
                             # theta0 = seq(-5,  3,length.out = 6), #theta0 = seq(0.2, 1.5,length.out = 3),
                             # thetamax = seq(-1, 5,length.out = 6), #thetamax = seq(1.6, 2.5,length.out = 3))
  }
  init_grid <- init_grid[init_grid$vmin < init_grid$vmax,]

  # Remove columns for fixed parameters
  init_grid <- init_grid[setdiff(names(init_grid), names(fixed))]
  # thetamax <- NULL # to omit a note because of an unbound variable
  # theta0 <- NULL # to omit a note because of an unbound variable
  # init_grid <- subset(init_grid, thetamax > theta0)
  init_grid <- unique(init_grid)

  ### If no grid-search is desired use mean of possible parameters #
  if (!grid_search) {
    init_grid <- as.data.frame(as.list(colMeans(init_grid)))
  }

  ##### span drifts with quadratic distance
  ###  We assume a different V (mean drift rate) for the different conditions --> nConds parameters
  if (nConds==1) {
    init_grid$v1 <- (init_grid$vmin+init_grid$vmax)/2
  } else {
    for (i in 0:(nConds-1)){
      init_grid[paste("v", i+1, sep="")] <- init_grid$vmin+(i/(nConds-1))^2*(init_grid$vmax-init_grid$vmin)
    }
  }
  init_grid$thetaMin <- init_grid$thetaMin - init_grid$t0 - init_grid$st0
  init_grid$thetaMax <- init_grid$thetaMax - init_grid$t0
  init_grid <- init_grid[init_grid$thetaMin<init_grid$thetaMax,]
  init_grid <- init_grid[init_grid$thetaMin>0,]

  if (sym_thetas) {
    init_grid[paste0("theta",(nRatings-1))] <- init_grid$thetaMin
    if (nRatings > 2) {
      for (i in (nRatings-2):1) {
        init_grid[paste("dtheta", i, sep="")] <- (init_grid$thetaMax -init_grid$thetaMin) / (nRatings-2)
      }
      cols_theta <- c(paste0("theta",(nRatings-1)), paste("dtheta", (nRatings-2):1, sep=""))
    } else {
      cols_theta <- c("theta1")
    }
  } else {
    init_grid[paste0(c("thetaUpper", "thetaLower"),(nRatings-1))] <- init_grid$thetaMin
    if (nRatings > 2) {
      for (i in (nRatings-2):1) {
        init_grid[paste(c("dthetaUpper", "dthetaLower"), i, sep="")] <-
          (init_grid$thetaMax -init_grid$thetaMin) / (nRatings-2)
      }
      cols_theta <- c(paste0(c("thetaLower"),(nRatings-1)), paste("dthetaLower", (nRatings-2):1, sep=""),
                      paste0(c("thetaUpper"),(nRatings-1)), paste("dthetaUpper", (nRatings-2):1, sep=""))
    } else {
      cols_theta <- c("thetaLower1", "thetaUpper1")
    }
  }
  parnames <- c('a', 'z', 'sz', paste("v", 1:nConds, sep=""),
                 'st0', 'sv', 't0', cols_theta)
  inits <- init_grid[, setdiff(parnames, names(fixed))]


  # remove init_grid
  rm(init_grid)

  ## Intermezzo: Setup cluster for parallelization   ####
  if (useparallel) {
    if (is.null(n.cores)) {
      n.cores <- detectCores()-1
    }
    cl <- makeCluster(type="SOCK", n.cores)
    clusterExport(cl, c("df", "fixed",
                        "nConds","nRatings", "sym_thetas", "precision"), envir = environment())
  }



  ### 2. Search initial grid before optimization  ####
  if (grid_search) {
    if (logging==TRUE) {
      logger::log_info(paste(length(inits[,1]), "...parameter sets to check"))
      logger::log_info(paste("data was compressed due to rounding rts by a factor of " , round(sum(df$n)/nrow(df), digits=2), "."))
      logger::log_info(paste("fitting data got ", nrow(df), " rows"))

      t00 <- Sys.time()
      logger::log_info("Searching initial values ...")
    }

    if (useparallel) {
      logL <-
        parApply(cl, inits, MARGIN=1,
                 function(p) try(neglikelihood_DDMConf_bounded(p, df,nConds, nRatings, fixed,sym_thetas, precision, st0stepsize),
                                 silent=TRUE))
      #stopCluster(cl)
    } else {
      logL <-
        apply(inits, MARGIN = 1,
              function(p) try(neglikelihood_DDMConf_bounded(p, df, nConds, nRatings, fixed,sym_thetas, precision, st0stepsize),
                              silent=TRUE))
    }
    logL <- as.numeric(logL)
    inits <- inits[order(logL),]
    if (logging) {
      logger::log_success(paste("Initial grid search took...",as.character(round(as.double(difftime(Sys.time(),t00,units = "mins")), 2))," mins"))
    }
  } else {
    logL <- NULL
  }

                    # a,  z, sz,  v1, v2,....,,   st0,             sv, t0,  thetaLower_nRatings-1, dthetaLower2.., thetaUpper1... (or theta1,...)
  lower_optbound <- c(0,  0,  0,  rep(0, nConds), minst0,          0,  0,   rep(c(0,     rep(0, nRatings-2)), 2-as.numeric(sym_thetas)))[!(parnames %in% names(fixed))]
  upper_optbound <- c(Inf,1,  1,  rep(Inf,nConds),max(max(df$rt)+1, 8),Inf,Inf, rep(Inf, (2-as.numeric(sym_thetas))*(nRatings-1))          )[!(parnames %in% names(fixed))]



  #### 3. Optimization ####
  if (logging==TRUE) {
    logger::log_info("Start fitting ... ")
  }
  if (!useparallel || (opts$nAttempts==1)) {
    noFitYet <- TRUE
    for (i in 1:opts$nAttempts){
      start <- c(t(inits[i,]))
      names(start) <- names(inits)
      for (l in 1:opts$nRestarts){
        try(m <- bobyqa(par = start,
                        fn = neglikelihood_DDMConf_bounded,
                        lower = lower_optbound, upper = upper_optbound,
                        data=df, nConds=nConds, nRatings=nRatings,
                        fixed = fixed, sym_thetas=sym_thetas, precision=precision,
                        st0stepsize=st0stepsize,
                        control = list(maxfun=opts$maxfun,
                                       rhobeg = min(0.02, 0.2*max(abs(start))),
                                       npt = length(start)+5)))
        ## rhobeg should be: about 0.1*(greatest expected change in parameters --> <= 1-2 (for a, thetas or v's) )
        ##                   smaller than min(abs(upper-lower)) = min(1, restr_tau)
        ##                   --> so we use min(0.2*restr_tau, 0.2, 0.2*max(abs(par))).
        ##                   Default would be: min(0.95, 0.2*max(abs(par))), respectively 0.2*max(upper_optbound-lower_optbound)
        ## rhoend: use default of 1e-6*rhobeg
        if (exists("m") && !inherits(m, "try-error")){
          m$value <- m$fval
        }

        if (logging==TRUE) {
          logger::log_info(paste("Finished attempt No.", i, " restart no. ", l))
        }
        if (!exists("m") || inherits(m, "try-error")){
          if (logging==TRUE) {
            logger::log_error(paste("No fit obtained at attempt No.", i))
            logger::log_error(paste("Used parameter set", paste(start, sep="", collapse=" "), sep=" ", collapse = ""))
          }
          break
        }
        if (exists("m") && is.list(m)){
          if (noFitYet) {
            fit <- m
            noFitYet <- FALSE
            if (logging==TRUE) {
              logger::log_info(paste("First fit obtained at attempt No.", i))
              attempt <- i
              save(logL, inits,  df,fit, attempt,file=filename)
            }
            start <- fit$par
          } else if (m$value < fit$value) {
            fit <- m
            if (logging==TRUE) {
              logger::log_info(paste("New fit at attempt No.", i, " restart no. ", l))
              attempt <- i
              save(logL, inits,  df,fit, attempt,file=filename)
            }
            start <- fit$par
          } # end of if better value
        }   # end of if we got a optim-result at all
      }     # end of for restarts
    }       # end of for initial start values
  } else {  # if useparallel
    starts <- inits[(1:opts$nAttempts),]
    parnames <- names(starts)

    optim_node <- function(start) { # define optim-routine to run on each node
      noFitYet <- TRUE
      start <- c(t(start))
      names(start) <- parnames
      for (l in 1:opts$nRestarts){
        try(m <- bobyqa(par = start,
                        fn = neglikelihood_DDMConf_bounded,
                        lower = lower_optbound, upper = upper_optbound,
                        data=df,  nConds=nConds, nRatings=nRatings,
                        fixed = fixed, sym_thetas=sym_thetas, precision=precision,
                        st0stepsize = st0stepsize,
                        control = list(maxfun=opts$maxfun,
                                       rhobeg = min(0.02, 0.2*max(abs(start))),
                                       npt = length(start)+5)),
            silent=TRUE)
        ## rhobeg should be: about 0.1*(greatest expected change in parameters --> <= 1-2 (for a, thetas or v's) )
        ##                   smaller than min(abs(upper-lower)) = min(1, restr_tau)
        ##                   --> so we use min(0.2*restr_tau, 0.2, 0.2*max(abs(par))).
        ##                   Default would be: min(0.95, 0.2*max(abs(par))), respectively 0.2*max(upper_optbound-lower_optbound)
        ## rhoend: use default of 1e-6*rhobeg
        if (exists("m") && !inherits(m, "try-error")){
          m$value <- m$fval
        }
        if (!exists("m") || inherits(m, "try-error")){
          break
        }
        if (exists("m") && is.list(m)){
          if (noFitYet) {
            fit <- m
            noFitYet <- FALSE
            start <- fit$par
          } else if (m$value < fit$value) {
            fit <- m
            start <- fit$par
          }
        }
      }
      if (exists("fit") && is.list(fit)){
        return(c(fit$value,fit$par))
      } else {
        # If optimization broke return starting values and start-negloglik
        return(c(start, neglikelihood_DDMConf_bounded(start, df,nConds, nRatings, fixed,sym_thetas, precision, st0stepsize)))
      } # end of node-function
    }
    clusterExport(cl, c("parnames", "opts", "optim_node" ), envir = environment())
    clusterExport(cl, c("lower_optbound", "upper_optbound"), envir = environment())

    optim_outs <- parApply(cl, starts,MARGIN=1, optim_node )
    stopCluster(cl)
    optim_outs <- t(optim_outs)
    best_res <- optim_outs[order(optim_outs[,1]),][1,]
    fit <- list(par = best_res[-1], value=best_res[1])
  }   # end of if-else useparallel

  #### 4. Wrap up results ####
  res <-  data.frame(matrix(nrow=1, ncol=0))
  if(exists("fit") && is.list(fit)){
    k <- length(fit$par)
    N <- nrow(df)
    p <- fit$par


    p <- c(t(p))
    names(p) <- names(inits)
    if (nRatings>2) {
      if (sym_thetas) {
        p[paste("theta", (nRatings-2):1, sep="")] <- c(t(p[paste("theta", (nRatings-1), sep="")])) + cumsum(c(t(p[grep(names(p), pattern = "dtheta", value=TRUE)])))
      } else {
        p[paste("thetaUpper", (nRatings-2):1, sep="")] <- c(t(p[paste("thetaUpper", (nRatings-1), sep="")])) + cumsum(c(t(p[grep(names(p), pattern = "dthetaUpper", value=TRUE)])))
        p[paste("thetaLower", (nRatings-2):1, sep="")] <- c(t(p[paste("thetaLower", (nRatings-1), sep="")])) + cumsum(c(t(p[grep(names(p), pattern = "dthetaLower", value=TRUE)])))
      }
      p <- p[ -grep(names(p), pattern="dtheta")]
    }
    if (!("sz" %in% names(fixed))) {
      if (!("z" %in% names(fixed))) {
        p['sz'] <- (min(p['z'], (1-p['z']))*2)*p["sz"] ##
      } else {
        p['sz'] <- (min(fixed[['z']], (1-fixed[['z']]))*2)*p["sz"] ##
      }
    }
    res <-   data.frame(matrix(nrow=1, ncol=length(p)))
    res[1,] <- p
    names(res) <- names(p)      # a,  z, sz,v1, v2,....,,   st0, sv, t0, thetaLower1,dthetaLower2-4,   thetaUpper1,dthetaUpper2-4,    tau


    if (!is.null(used_cats)) {
      # If some rating categories are not used, we fit less thresholds numerically and fill up the
      # rest by the obvious best-fitting thresholds (e.g. +/- Inf for the lowest/highest...)
      res <- fill_thresholdsDDM(res, used_cats, actual_nRatings)
      nRatings <- actual_nRatings
      k <- ncol(res)
    }

    if (sym_thetas) {
      parnames <- c(paste("v", 1:nConds, sep=""), 'sv', 'a', 'z', 'sz', 't0','st0', paste("theta", 1:(nRatings-1), sep=""))
    } else {
      parnames <- c(paste("v", 1:nConds, sep=""), 'sv', 'a', 'z', 'sz', 't0','st0', paste("thetaLower", 1:(nRatings-1), sep=""), paste("thetaUpper", 1:(nRatings-1), sep=""))
    }
    res <- res[,setdiff(parnames, names(fixed))]
    if (length(fixed)>=1) {
      res <- cbind(res, as.data.frame(fixed))
    }
    res$fixed <- paste(c("sym_thetas", names(fixed)), c(sym_thetas,fixed), sep="=", collapse = ", ")
    res$negLogLik <- fit$value
    res$N <- N
    res$k <- k
    res$BIC <-  2 * fit$value + k * log(N)
    res$AICc <- 2 * fit$value + k * 2 + 2*k*(k-1)/(N-k-1)
    res$AIC <- 2 * fit$value + k * 2
    if (logging==TRUE) {
      logger::log_success("Done fitting and autosaved results")
      save(logL, df, res, file=filename)
    }
  }
  return(res)
}




neglikelihood_DDMConf_bounded <-   function(p, data, nConds, nRatings, fixed, sym_thetas=FALSE, precision=1e-5, st0stepsize=0.001)
{
  # get parameter vector back from real transformations
  paramDf <-   data.frame(matrix(nrow=1, ncol=length(p)))
  paramDf[1,] <- p
  names(paramDf) <- names(p)
  if (nRatings > 2) {
    if (sym_thetas) {
      paramDf[paste("theta", (nRatings-2):1, sep="")] <- c(t(paramDf[paste("theta", (nRatings-1), sep="")])) + cumsum(c(t(paramDf[grep(names(paramDf), pattern = "dtheta", value=TRUE)])))
    } else {
      paramDf[paste("thetaUpper", (nRatings-2):1, sep="")] <- c(t(paramDf[paste("thetaUpper", (nRatings-1), sep="")])) + cumsum(c(t(paramDf[grep(names(paramDf), pattern = "dthetaUpper", value=TRUE)])))
      paramDf[paste("thetaLower", (nRatings-2):1, sep="")] <- c(t(paramDf[paste("thetaLower", (nRatings-1), sep="")])) + cumsum(c(t(paramDf[grep(names(paramDf), pattern = "dthetaLower", value=TRUE)])))
    }
    paramDf <- paramDf[ -grep(names(paramDf), pattern="dtheta")]
  }
  if (length(fixed)>=1) {
    paramDf <- cbind(paramDf, as.data.frame(fixed))
  }
  paramDf['sz'] <- (min(paramDf['z'], 1-paramDf['z'])*2)*paramDf['sz']
  #paramDf[grep(names(paramDf), pattern = "thetaLower", value=TRUE)] <- paramDf$a - paramDf[grep(names(paramDf), pattern = "thetaLower", value=TRUE)]



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
  vary_s <-   length(grep(pattern = "^s[0-9]", names(paramDf), value = T))>1
  if (vary_s){
    S <- c(t((paramDf[,paste("s",1:(nConds), sep = "")])))
  } else {
    if ("s" %in% names(paramDf)) {
      S <- rep(paramDf$s, nConds)
    } else {
      S <- rep(1, nConds)
    }
  }
  ## Recover confidence thresholds
  if (sym_thetas) {
    thetas_upper <- c(.Machine$double.eps, t(paramDf[,paste("theta",(nRatings-1):1, sep = "")]), 1e+64)
    thetas_lower <- c(.Machine$double.eps, t(paramDf[,paste("theta",(nRatings-1):1, sep = "")]), 1e+64)
  } else {
    thetas_upper <- c(.Machine$double.eps, t(paramDf[,paste("thetaUpper",(nRatings-1):1, sep = "")]), 1e+64)
    thetas_lower <- c(.Machine$double.eps, t(paramDf[,paste("thetaLower",(nRatings-1):1, sep="")]), 1e+64)
  }

  ## Compute the row-wise likelihood of observations
  data <-data %>% mutate(th1 = case_when(.data$response ==  1 ~ thetas_upper[(nRatings + 1 - .data$rating)],
                                          .data$response== -1 ~ thetas_lower[(nRatings + 1 - .data$rating)]),
                         th2 = case_when(.data$response ==  1 ~ thetas_upper[(nRatings + 2 - .data$rating)],
                                          .data$response== -1 ~ thetas_lower[(nRatings + 2 - .data$rating)]),
                         M_drift = V[.data$condition]*.data$stimulus,
                         SV = SV[.data$condition],
                         S = S[.data$condition])
  probs <- with(data, dDDMConf(rt,  response, th1, th2,
                                  a=paramDf$a,
                                  v = M_drift,
                                  t0 = paramDf$t0, z = paramDf$z, sz = paramDf$sz, st0=paramDf$st0,
                                  sv = SV, s = S,
                                  z_absolute = FALSE,
                                  precision = precision, stop_on_error = TRUE,
                                  stop_on_zero=FALSE, st0stepsize = st0stepsize))

  ## Produce output as log-Likelihood
  if (any(is.na(probs))) return(1e12)
  # if (any(probs<=0)) {
  #   return(1e12)
  # }
  probs[probs==0] <- 1e-322
  if ("n" %in% names(data)) {
    negloglik <- -sum(log(probs)*data$n)
  } else {
    negloglik <- -sum(log(probs))
  }
  return(negloglik)
}

