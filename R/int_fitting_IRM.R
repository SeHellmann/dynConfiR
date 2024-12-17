fittingIRM <- function(df, nConds, nRatings, fixed, sym_thetas, time_scaled,
                          grid_search, init_grid=NULL, optim_method, opts,
                          logging, filename,
                          useparallel, n.cores,
                          used_cats, actual_nRatings, precision){
  ## Be sure that the parallel cluster is stopped if anything happens (error or user interupt)
  on.exit(try(stopCluster(cl), silent = TRUE))

  fitted_weights <- NULL
  if (time_scaled) {   ## Incorporate fixed weight parameters
    max_par_weight <- 1
    fixed_ws <- grep(pattern="w", names(fixed), value=TRUE)
    if (length(fixed_ws)>1) { ## If 2 are fixed the third weight is also fixed, since they have to sum to 1
      if ((length(fixed_ws)==3) && (sum(as.numeric(fixed[fixed_ws]))!= 1)) stop("Provided weights do not sum to 1!")
      if ((length(fixed_ws)==2)) fixed[[setdiff(c("wx","wrt", "wint"),fixed_ws)]] <- 1- sum(as.numeric(fixed[fixed_ws]))
      fitted_weights <- NULL
      wx <- 0
      wrt <- 0
    }
    if (length(fixed_ws)==1) { ## If only 1 is supplied, one may be fit
      max_par_weight <- fixed[[fixed_ws]]
      fitted_weights <- setdiff(c("wx","wrt", "wint"),fixed_ws)[1]
      wx <- 0
      wrt <- 0
      assign(fitted_weights, c(0.1, 0.5, 0.9)) # Range of fitted parameter: (0,1),
      # but this will be multiplied with max_fit_weight after optimisation

    }
    if (length(fixed_ws)==0) { ## If no weight is supplied, two have to be fit
      wx <- c(0.025, 0.3, 0.8)
      ## To keep the optimization box-constraint, we fit the second parameter
      ## as proportion of the rest of 1 after considering wx, i.e.
      ## wrt(true weight) = (1-wx)*wrt(fitted parameter)
      wrt <- c(0.05, 0.5, 0.9)
      fitted_weights <- c("wx","wrt")
    }
  }
  ### 1. Generate initial grid for grid search over possible parameter sets ####
  #### Create grid ####
  mint0 <-  min(df$rt)
  if (is.null(init_grid)) {
    #if (mint0 < 0.2) {mint0 <- 0.2}
      if (time_scaled) {
        init_grid <- expand.grid(vmin = c(0.01, 0.1,0.8),         ### vmin = drift rate in first condition \in (0,\infty)]
                                 vmax = c(1, 2, 3.8, 5),          ### vmax = mean drift rate in last condition \in (\vmin,\infty)]
                                 a = c(0.5, 1.2, 2.8),
                                 b = c(0.5, 1.2, 2.8),
                                 theta0 = c(0.3,1.2, 1.7, 2.5),  ### theta0 = lowest threshold for confidence rating (in difference from threshold / a)
                                 thetamax = c(1, 3, 3, 4),       ### thetamax = highest threshold for confidence rating (in distance from threshold / a)
                                 t0 = seq(max(mint0-0.2, 0), max(mint0-0.1, 0)), ### t0 = minimal non-decision time
                                 st0 = seq(0.07, 1, length.out=3),    ### st0 = range of (uniform dist) non-decision time
                                 wx = wx,      ### coeff for BoE in conf= wx *(b-xj) + (wrt* 1//sqrt(t)) + (wint* (b-xj)/sqrt(t))
                                 wrt = wrt)      ### coeff for time in conf= wx *(b-xj) + (wrt* 1//sqrt(t)) + (wint* (b-xj)/sqrt(t))
      } else {
        init_grid <- expand.grid(vmin = c(0.01, 0.1,0.8),         ### vmin = drift rate in first condition \in (0,\infty)]
                                 vmax = c(1, 2, 3.8, 5),          ### vmax = mean drift rate in last condition \in (\vmin,\infty)]
                                 a = c(0.5, 1.2, 2.8),
                                 b = c(0.5, 1.2, 2.8),
                                 theta0 = c(0.3,1.2, 1.7, 2.5),  ### theta0 = lowest threshold for confidence rating (in difference from threshold / a)
                                 thetamax = c(1, 3, 3, 4),       ### thetamax = highest threshold for confidence rating (in distance from threshold / a)
                                 t0 = seq(max(mint0-0.2, 0), max(mint0-0.1, 0)), ### t0 = minimal non-decision time
                                 st0 = seq(0.07, 1, length.out=3))    ### st0 = range of (uniform dist) non-decision time
      }
    }
  init_grid <- init_grid[init_grid$theta0 < init_grid$thetamax,]
  # Remove column for weight that will not be fitted
  init_grid <- init_grid[c(setdiff(names(init_grid), c("wx","wrt")), fitted_weights)]
  # Remove columns for fixed parameters
  init_grid <- init_grid[setdiff(names(init_grid), names(fixed))]
  init_grid <- unique(init_grid)

  ### If no grid-search is desired use mean of possible parameters #
  if (!grid_search) {
    init_grid <- as.data.frame(as.list(colMeans(init_grid)))
  }

  #### 1.1. For Nelder-Mead transform all parameters to real values ####
  if (optim_method=="Nelder-Mead") {
    ## change parametrisation (should be on the whole real line) and
    #  span V-parameters between vmin and vmax equidistantly for all conditions
    inits <- data.frame(matrix(data=NA, nrow= nrow(init_grid),
                               ncol = nConds))
    for (i in 0:(nConds-1)){
      if (nConds == 1) {
        inits[,1] <- log((init_grid$vmin+init_grid$vmax)/2)
      } else {
        inits[,i+1] <- log(init_grid$vmin+(i/(nConds-1))^3*(init_grid$vmax-init_grid$vmin)) ###  We assume a different V (mean drift rate) for the different conditions --> nConds parameters
      }
    }
    if (!("a" %in% names(fixed))) inits <- cbind(inits, log(init_grid$a))
    if (!("b" %in% names(fixed))) inits <- cbind(inits, log(init_grid$b))
    if (!("t0" %in% names(fixed))) inits <- cbind(inits, log(init_grid$t0))
    if (!("st0" %in% names(fixed))) inits <- cbind(inits, log(init_grid$st0))
    if (time_scaled) {
      for (par in fitted_weights) {
        inits <- cbind(inits, qnorm(init_grid[[par]]))
      }
    }

    inits <- cbind(inits, init_grid$theta0)
    if (nRatings > 2) {
      for (i in 1:(nRatings-2)) {
        inits <- cbind(inits, log((init_grid$thetamax-init_grid$theta0)/(nRatings-2)))
      }
    }
    if (!sym_thetas) {
      inits <- cbind(inits, init_grid$theta0)
      if (nRatings > 2) {
        for (i in 1:(nRatings-2)) {
          inits <- cbind(inits, log((init_grid$thetamax-init_grid$theta0)/(nRatings-2)))
        }
      }
      cols_theta <- c('thetaLower1', rep(paste("dthetaLower", 2:(nRatings-1), sep=""), times=nRatings>2),
                      'thetaUpper1', rep(paste("dthetaUpper", 2:(nRatings-1), sep=""), times=nRatings>2))
    } else {
      cols_theta <- c("theta1", rep(paste("dtheta", 2:(nRatings-1), sep=""), times=nRatings>2))
    }

    ## replace all +-Inf with big/tiny numbers
    inits[inits==Inf]<- 1e6
    inits[inits==-Inf]<- -1e6
    parnames <- c(paste("v", 1:nConds, sep=""), 'a', 'b', 't0', 'st0', fitted_weights, cols_theta)
    names(inits) <- setdiff(parnames, names(fixed))

  } else {

    ##### 1.2. For box-constraint optimisation algorithm span drifts and thresholds equidistantly
    if (nConds==1) {
      init_grid$v1 <- (init_grid$vmin+init_grid$vmax)/2
    } else {
      for (i in 0:(nConds-1)){
        init_grid[paste("v", i+1, sep="")] <- init_grid$vmin+(i/(nConds-1))^3*(init_grid$vmax-init_grid$vmin)
      }
    }

    if (sym_thetas) {
      init_grid["theta1"] <- init_grid$theta0
      if (nRatings > 2) {
        for (i in 2:(nRatings-1)) {
          init_grid[paste("dtheta", i, sep="")] <- (init_grid$thetamax-init_grid$theta0)/(nRatings-2)
        }
        cols_theta <- c("theta1", paste("dtheta", 2:(nRatings-1), sep=""))
      } else {
        cols_theta <- c("theta1")
      }
    } else {
      init_grid[c("thetaUpper1", "thetaLower1")] <- init_grid$theta0
      if (nRatings > 2) {
        for (i in 2:(nRatings-1)) {
          init_grid[paste(c("dthetaUpper", "dthetaLower"), i, sep="")] <- (init_grid$thetamax-init_grid$theta0)/(nRatings-2)
        }
        cols_theta <- c('thetaLower1', paste("dthetaLower", 2:(nRatings-1), sep=""),
                        'thetaUpper1', paste("dthetaUpper", 2:(nRatings-1), sep=""))
      } else {
        cols_theta <- c("thetaLower1", "thetaUpper1")
      }
    }
    parnames <- c('a', 'b', 't0', 'st0', fitted_weights, paste("v", 1:nConds, sep=""), cols_theta)
    inits <- init_grid[, setdiff(parnames, names(fixed))]
  }
  # remove init_grid
  rm(init_grid)


  ## Intermezzo: Setup cluster for parallelization   ####
  if (useparallel) {
    if (is.null(n.cores)) {
      n.cores <- detectCores()-1
    }
    cl <- makeCluster(type="SOCK", n.cores)
    clusterExport(cl, c("df", "time_scaled",
                        "nConds","nRatings", "sym_thetas", "fixed", "fitted_weights"), envir = environment())
  }



  ### 2. Search initial grid before optimization  ####
  if (grid_search) {
    if (logging==TRUE) {
      logger::log_info(paste(length(inits[,1]), "...parameter sets to check"))
      logger::log_info(paste("data got ", nrow(df), " rows"))
      t00 <- Sys.time()
      logger::log_info("Searching initial values ...")
    }

    if (optim_method =="Nelder-Mead") {
      if (useparallel) {
        logL <-
          parApply(cl, inits, MARGIN=1,
                   function(p) try(neglikelihood_IRM_free(p, df,
                                                             time_scaled, nConds, nRatings, fixed, fitted_weights, sym_thetas, precision),
                                   silent=TRUE))
        #stopCluster(cl)
      } else {
        logL <-
          apply(inits, MARGIN = 1,
                function(p) try(neglikelihood_IRM_free(p, df, time_scaled, nConds, nRatings, fixed, fitted_weights, sym_thetas, precision),
                                silent = TRUE))
      }
    } else {
      if (useparallel) {
        logL <-
          parApply(cl, inits, MARGIN=1,
                   function(p) try(neglikelihood_IRM_bounded(p, df, time_scaled, nConds, nRatings,fixed, fitted_weights,  sym_thetas, precision),
                                   silent=TRUE))
        #stopCluster(cl)
      } else {
        logL <-
          apply(inits, MARGIN = 1,
                function(p) try(neglikelihood_IRM_bounded(p, df, time_scaled, nConds, nRatings, fixed, fitted_weights, sym_thetas, precision),
                                silent=TRUE))
      }
    }
    logL <- as.numeric(logL)
    inits <- inits[order(logL),]
    if (logging==TRUE) {
      logger::log_success(paste("Initial grid search took...",as.character(round(as.double(difftime(Sys.time(),t00,units = "mins")), 2))," mins"))
    }
  } else {
    logL <- NULL
  }

  if (optim_method!="Nelder-Mead") {
                      # a,  b,  t0, st0, wrt and/or wint, (if needed (time_scaled=TRUE)),      v1, v2,....,,   thetaLower1, dthetaLower2.., thetaUpper1... (or theta1,...)
    lower_optbound <- c(0,  0,  0,  0,   rep(0,length(fitted_weights)*as.numeric(time_scaled)),rep(0, nConds), rep(rep(0, nRatings-1), 2-as.numeric(sym_thetas)))[!(parnames %in% names(fixed))]
    upper_optbound <- c(Inf,Inf,Inf,Inf, rep(1,length(fitted_weights)*as.numeric(time_scaled)),rep(Inf,nConds),rep(Inf, (2-as.numeric(sym_thetas))*(nRatings-1)))[!(parnames %in% names(fixed))]
  }


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
        if (optim_method == "Nelder-Mead") {
          try(m <- optim(par = start,
                         fn = neglikelihood_IRM_free,
                         data=df, time_scaled=time_scaled, nConds=nConds, nRatings=nRatings,
                         fixed=fixed, fitted_weights=fitted_weights,
                         sym_thetas=sym_thetas, precision=precision,
                         method="Nelder-Mead",
                         control = list(maxit = opts$maxit, reltol = opts$reltol)))
        } else if (optim_method =="bobyqa") {
          start <- pmax(pmin(start, upper_optbound-1e-6), lower_optbound+1e-6)
          try(m <- bobyqa(par = start,
                          fn = neglikelihood_IRM_bounded,
                          lower = lower_optbound, upper = upper_optbound,
                          data=df, time_scaled=time_scaled, nConds=nConds, nRatings=nRatings,
                          fixed=fixed, fitted_weights=fitted_weights,
                          sym_thetas=sym_thetas, precision=precision,
                          control = list(maxfun=opts$maxfun,
                                         rhobeg = min(0.2, 0.2*max(abs(start))),
                                         npt = length(start)+5)))
          ## rhobeg should be: about 0.1*(greatest expected change in parameters --> <= 1-2 (for a, thetas or v's) )
          ##                   Default would be: min(0.95, 0.2*max(abs(par)))
          ## rhoend: use default of 1e-6*rhobeg
          if (exists("m") && !inherits(m, "try-error")){
            m$value <- m$fval
          }
        } else if (optim_method=="L-BFGS-B") {  ### ToDo: use dfoptim or pracma::grad as gradient!
          start <- pmax(pmin(start, upper_optbound-1e-6), lower_optbound+1e-6)
          try(m <- optim(par = start,
                         fn = neglikelihood_IRM_bounded,
                         lower = lower_optbound, upper = upper_optbound,
                         data=df,  time_scaled=time_scaled, nConds=nConds, nRatings=nRatings,
                         fixed=fixed, fitted_weights=fitted_weights,
                         sym_thetas=sym_thetas, precision=precision,
                         method="L-BFGS-B",
                         control = list(maxit = opts$maxit, factr = opts$factr)))
        } else {
          stop(paste("Not implemented or unknown method: ", optim_method, ". Use 'bobyqa', Nelder-Mead' or 'L-BFGS-B' instead.", sep=""))
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
            names(start) <- names(inits)
          } else if (m$value < fit$value) {
            fit <- m
            if (logging==TRUE) {
              logger::log_info(paste("New fit at attempt No.", i, " restart no. ", l))
              attempt <- i
              save(logL, inits,  df,fit, attempt,file=filename)
            }
            start <- fit$par
            names(start) <- names(inits)
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
        if (optim_method == "Nelder-Mead") {
          try(m <- optim(par = start,
                         fn = neglikelihood_IRM_free,
                         data=df,  time_scaled=time_scaled, nConds=nConds, nRatings=nRatings,
                         fixed=fixed, fitted_weights=fitted_weights,
                         sym_thetas=sym_thetas, precision=precision,
                         method="Nelder-Mead",
                         control = list(maxit = opts$maxit, reltol = opts$reltol)))
        } else if (optim_method =="bobyqa") {
          start <- pmax(pmin(start, upper_optbound-1e-6), lower_optbound+1e-6)
          try(m <- bobyqa(par = start,
                          fn = neglikelihood_IRM_bounded,
                          lower = lower_optbound, upper = upper_optbound,
                          data=df,  time_scaled=time_scaled, nConds=nConds, nRatings=nRatings,
                          fixed=fixed, fitted_weights=fitted_weights,
                          sym_thetas=sym_thetas, precision=precision,
                          control = list(maxfun=opts$maxfun,
                                         rhobeg = min(0.2, 0.2*max(abs(start))),
                                         npt = length(start)+5)),
              silent=TRUE)
          ## rhobeg should be: about 0.1*(greatest expected change in parameters --> <= 1-2 (for a, thetas or v's) )
          ##                   Default would be: min(0.95, 0.2*max(abs(par)))
          ## rhoend: use default of 1e-6*rhobeg
          if (exists("m") && !inherits(m, "try-error")){
            m$value <- m$fval
          }
        } else if (optim_method=="L-BFGS-B") {  ### ToDo: use dfoptim or pracma::grad as gradient!
          start <- pmax(pmin(start, upper_optbound-1e-6), lower_optbound+1e-6)
          try(m <- optim(par = start,
                         fn = neglikelihood_IRM_bounded,
                         lower = lower_optbound, upper = upper_optbound,
                         data=df,  time_scaled=time_scaled, nConds=nConds, nRatings=nRatings,
                         fixed=fixed, fitted_weights=fitted_weights,
                         sym_thetas=sym_thetas, precision=precision,
                         method="L-BFGS-B",
                         control = list(maxit = opts$maxit, factr = opts$factr)))
        } else {
          stop(paste("Not implemented or unknown method: ", optim_method, ". Use 'bobyqa', Nelder-Mead' or 'L-BFGS-B' instead.", sep=""))
        }
        if (!exists("m") || inherits(m, "try-error")){
          break
        }
        if (exists("m") && is.list(m)){
          if (noFitYet) {
            fit <- m
            noFitYet <- FALSE
            start <- fit$par
            names(start) <- parnames
          } else if (m$value < fit$value) {
            fit <- m
            start <- fit$par
            names(start) <- parnames
          }
        }
      }
      if (exists("fit") && is.list(fit)){
        return(c(fit$value,fit$par))
      } else {
        return(rep(NA, length(start)+1))
      } # end of node-function
    }
    clusterExport(cl, c("parnames", "opts", "optim_method","optim_node" ), envir = environment())
    if (optim_method!="Nelder-Mead") {
      clusterExport(cl, c("lower_optbound", "upper_optbound"), envir = environment())
    }
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
    N <- sum(df$n)

    p <- c(t(fit$par))
    if (optim_method=="Nelder-Mead") {
      names(p) <- names(inits)
      res[,paste("v",1:(nConds), sep="")] <- exp(p[1:(nConds)])
      if (!("a" %in% names(fixed))) res$a <- exp(p[["a"]])
      if (!("b" %in% names(fixed))) res$b <- exp(p[["b"]])
      if (!("t0" %in% names(fixed))) res$t0 <- exp(p[["t0"]])
      if (!("st0" %in% names(fixed))) res$st0 <- exp(p[["st0"]])

      if (time_scaled) {
        if (length(fitted_weights) == 1) {
          res[,fitted_weights] <- pnorm(p[[fitted_weights]])*(1-fixed[[grep("^w", names(fixed), value = TRUE)]])
          res[,grep("^w", names(fixed), value = TRUE)] <- fixed[[grep("^w", names(fixed), value = TRUE)]]
          res[, setdiff(c("wx", "wrt", "wint"),names(res))] <- 1- sum(as.numeric(res[grep("^w", names(res), value=TRUE)]))
        }
        if (length(fitted_weights)==2) {
          res$wx <- pnorm(p[["wx"]])
          res$wrt <- pnorm(p[["wrt"]])*(1-res$wx)
          res$wint <- 1 - res$wx - res$wrt
        }
      } else {
        res$wx <- 1
        res$wrt <- 0
        res$wint <- 0
      }
      if (nRatings > 2) {
        if (sym_thetas) {
          res[,paste("theta",1:(nRatings-1), sep="")] <- cumsum(c(p[["theta1"]], exp(p[paste0("dtheta", 2:(nRatings-1))])))
        } else {
          res[,paste("thetaUpper",1:(nRatings-1), sep="")] <- cumsum(c(p[["thetaUpper1"]], exp(p[paste0("dthetaUpper", 2:(nRatings-1))])))
          res[,paste("thetaLower",1:(nRatings-1), sep="")] <- cumsum(c(p[["thetaLower1"]], exp(p[paste0("dthetaLower", 2:(nRatings-1))])))
        }
      } else {
        if (sym_thetas) {
          res[,paste("theta",1:(nRatings-1), sep="")] <- p[["theta1"]]
        } else {
          res[,paste("thetaUpper",1:(nRatings-1), sep="")] <- p[["thetaUpper1"]]
          res[,paste("thetaLower",1:(nRatings-1), sep="")] <- p[["thetaLower1"]]
        }
      }
    } else {
      res <-   data.frame(matrix(nrow=1, ncol=length(p)))
      res[1,] <- p
      names(res) <- names(inits)      # a, b, wrt and wint (maybe), v1, v2,....,, thetaLower1,dthetaLower2-4,   thetaUpper1,dthetaUpper2-4,
      if (time_scaled) {
        if (length(fitted_weights) == 1) {
          res[,fitted_weights] <- res[,fitted_weights]*(1-fixed[[grep("^w", names(fixed), value = TRUE)]])
          res[,grep("^w", names(fixed), value = TRUE)] <- fixed[[grep("^w", names(fixed), value = TRUE)]]
          res[, setdiff(c("wx", "wrt", "wint"),names(res))] <- 1- sum(as.numeric(res[grep("^w", names(res), value=TRUE)]))
        }
        if (length(fitted_weights)==2) {
          res$wrt <- res[["wrt"]]*(1-res$wx)
          res$wint <- 1 - res$wx - res$wrt
        }
        if (length(fitted_weights)==0) {
          res$wx <- fixed[['wx']]
          res$wrt <- fixed[['wrt']]
          res$wint <- fixed[['wint']]
        }
      } else {
        res$wx <- 1
        res$wrt <- 0
        res$wint <- 0
      }
      if (nRatings>2) {
        if (sym_thetas) {
          res[paste("theta", 2:(nRatings-1), sep="")] <- c(t(res['theta1'])) + cumsum(c(t(res[grep(names(res), pattern = "dtheta", value=TRUE)])))
        } else {
          res[paste("thetaUpper", 2:(nRatings-1), sep="")] <- c(t(res['thetaUpper1'])) + cumsum(c(t(res[grep(names(res), pattern = "dthetaUpper", value=TRUE)])))
          res[paste("thetaLower", 2:(nRatings-1), sep="")] <- c(t(res['thetaLower1'])) + cumsum(c(t(res[grep(names(res), pattern = "dthetaLower", value=TRUE)])))
        }
        res <- res[ -grep(names(res), pattern="dtheta")]
      }
    }
    if (!is.null(used_cats)) {
      # If some rating categories are not used, we fit less thresholds numerically and fill up the
      # rest by the obvious best-fitting thresholds (e.g. +/- Inf for the lowest/highest...)
      res <- fill_thresholds(res, used_cats, actual_nRatings, 0)
      k <- k+(as.numeric(!sym_thetas)+1)*(actual_nRatings-nRatings)
      nRatings <- actual_nRatings
    }
    if (length(fixed)>=1) {
      res <- cbind(res, as.data.frame(fixed))
    }
    if (res$a == "b") res$a <- res$b
    if (res$b == "a") res$b <- res$a

    if (sym_thetas) {
      res <- res[,c('a','b', 't0', 'st0', paste("v", 1:nConds, sep=""), paste("theta", 1:(nRatings-1), sep=""), 'wx', 'wrt', 'wint')]
    } else {
      res <- res[,c('a','b', 't0', 'st0', paste("v", 1:nConds, sep=""), paste("thetaLower", 1:(nRatings-1), sep=""), paste("thetaUpper", 1:(nRatings-1), sep=""), 'wx', 'wrt', 'wint')]
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


neglikelihood_IRM_free <-   function(p, data,  time_scaled,
                                     nConds, nRatings,
                                     fixed, fitted_weights, sym_thetas, precision)
{
  # get parameter vector back from real transformations
  paramDf <-  data.frame(matrix(nrow=1, ncol=0))
  if (length(fixed)>=1) {
    paramDf <- cbind(paramDf, as.data.frame(fixed))
  }
  paramDf[,paste("v",1:(nConds), sep="")] <- exp(p[1:(nConds)])
  if (!("a" %in% names(fixed))) paramDf$a <- exp(p[["a"]])
  if (!("b" %in% names(fixed))) paramDf$b <- exp(p[["b"]])
  if (paramDf$a == "b") paramDf$a <- paramDf$b
  if (paramDf$b == "a") paramDf$b <- paramDf$a
  if (!("t0" %in% names(fixed))) paramDf$t0 <- exp(p[["t0"]])
  if (!("st0" %in% names(fixed))) paramDf$st0 <- exp(p[["st0"]])

  if (time_scaled) {
    if (length(fitted_weights) == 1) {
      paramDf[,fitted_weights] <- pnorm(p[[fitted_weights]])*(1-fixed[[grep("^w", names(fixed), value = TRUE)]])
      paramDf[,grep("^w", names(fixed), value = TRUE)] <- fixed[[grep("^w", names(fixed), value = TRUE)]]
      paramDf[, setdiff(c("wx", "wrt", "wint"),names(paramDf))] <- 1- sum(as.numeric(paramDf[grep("^w", names(paramDf), value=TRUE)]))
    }
    if (length(fitted_weights)==2) {
      paramDf$wx <- pnorm(p[["wx"]])
      paramDf$wrt <- pnorm(p[["wrt"]])*(1-paramDf$wx)
      paramDf$wint <- 1 - paramDf$wx - paramDf$wrt
    }
  } else {
    paramDf$wx <- 1
    paramDf$wrt <- 0
    paramDf$wint <- 0
  }
  if (nRatings > 2) {
    if (sym_thetas) {
      paramDf[,paste("theta",1:(nRatings-1), sep="")] <- cumsum(c(p[["theta1"]], exp(p[paste0("dtheta", 2:(nRatings-1))])))
    } else {
      paramDf[,paste("thetaUpper",1:(nRatings-1), sep="")] <- cumsum(c(p[["thetaUpper1"]], exp(p[paste0("dthetaUpper", 2:(nRatings-1))])))
      paramDf[,paste("thetaLower",1:(nRatings-1), sep="")] <- cumsum(c(p[["thetaLower1"]], exp(p[paste0("dthetaLower", 2:(nRatings-1))])))
    }
  } else {
    if (sym_thetas) {
      paramDf[,paste("theta",1:(nRatings-1), sep="")] <- p[["theta1"]]
    } else {
      paramDf[,paste("thetaUpper",1:(nRatings-1), sep="")] <- p[["thetaUpper1"]]
      paramDf[,paste("thetaLower",1:(nRatings-1), sep="")] <- p[["thetaLower1"]]
    }
  }
  if (any(is.infinite(t(paramDf))) || any(is.na(t(paramDf)))) {
    return(1e12)
  }
  negloglik <- -LogLikRM(data, paramDf, "IRM", time_scaled, precision=precision)

  return(negloglik)
}




neglikelihood_IRM_bounded <-   function(p, data, time_scaled,
                                        nConds, nRatings,
                                        fixed, fitted_weights,
                                        sym_thetas, precision)
{
  # get parameter vector back from real transformations
  paramDf <-   data.frame(matrix(nrow=1, ncol=length(p)))
  paramDf[1,] <- p
  names(paramDf) <- names(p)
  if (nRatings > 2) {
    if (sym_thetas) {
      paramDf[paste("theta", 2:(nRatings-1), sep="")] <- c(t(paramDf['theta1'])) + cumsum(c(t(paramDf[grep(names(paramDf), pattern = "dtheta", value=TRUE)])))
    } else {
      paramDf[paste("thetaUpper", 2:(nRatings-1), sep="")] <- c(t(paramDf['thetaUpper1'])) + cumsum(c(t(paramDf[grep(names(paramDf), pattern = "dthetaUpper", value=TRUE)])))
      paramDf[paste("thetaLower", 2:(nRatings-1), sep="")] <- c(t(paramDf['thetaLower1'])) + cumsum(c(t(paramDf[grep(names(paramDf), pattern = "dthetaLower", value=TRUE)])))
    }
    paramDf <- paramDf[ -grep(names(paramDf), pattern="dtheta")]
  }

  if (length(fixed)>=1) {
    paramDf <- cbind(paramDf, as.data.frame(fixed))
  }
  if (paramDf$a == "b") paramDf$a <- paramDf$b
  if (paramDf$b == "a") paramDf$b <- paramDf$a


  if (!time_scaled) {
    paramDf$wx <- 1
    paramDf$wrt <- 0
    paramDf$wint <- 0
  } else {
    if (length(fitted_weights) == 1) {
      paramDf[,fitted_weights] <- paramDf[,fitted_weights]*(1-fixed[[grep("^w", names(fixed), value = TRUE)]])
      paramDf[,grep("^w", names(fixed), value = TRUE)] <- fixed[[grep("^w", names(fixed), value = TRUE)]]
      paramDf[, setdiff(c("wx", "wrt", "wint"),names(paramDf))] <- 1- sum(as.numeric(paramDf[grep("^w", names(paramDf), value=TRUE)]))
    }
    if (length(fitted_weights)==2) {
      paramDf$wrt <- paramDf[["wrt"]]*(1-paramDf$wx)
      paramDf$wint <- 1 - paramDf$wx - paramDf$wrt
    }
  }
  negloglik <- -LogLikRM(data, paramDf, "IRM", time_scaled, precision=precision)
  return(negloglik)
}

