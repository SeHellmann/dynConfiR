fittingMTLNR <- function(df, nConds, nRatings, fixed, sym_thetas,
                       grid_search, init_grid=NULL, optim_method, opts,
                       logging, filename,
                       useparallel, n.cores,
                       used_cats, actual_nRatings, precision){
  ## Be sure that the parallel cluster is stopped if anything happens (error or user interupt)
  on.exit(try(stopCluster(cl), silent = TRUE))

  ### 1. Generate initial grid for grid search over possible parameter sets ####
  #### Create grid ####
  mint0 <-  min(df$rt)
  if (is.null(init_grid)) {
    # (mint0 < 0.2) {mint0 <- 0.2}
    init_grid <- expand.grid(vmin = c(0.1, 0.8, 3),         ### vmin = drift rate in first condition \in (0,\infty)]
                             vmax = c( 1, 3, 5),          ### vmax = mean drift rate in last condition \in (\vmin,\infty)]
                             mu_d1 = c(-0.5,0.5),
                             mu_d2 = c(-0.5, 0.5),
                             s_v1= c(0.5,3),
                             s_v2= c(0.5, 3),
                             s_d1= c(0.5, 3),
                             s_d2= c(0.5, 3),
                             rho_d= c(-0.4, 0.4),
                             rho_v= c(-0.4, 0.4),
                             t0 = c(0.2, 0.7), ### t0 = minimal non-decision time (proportional to minimum observed RT)
                             st0 = seq(0.1, 1.2, length.out=3))    ### st0 = range of (uniform dist) non-decision time
    # init_grid <- expand.grid(vmin = c(0.01, 0.1,0.8, 1.5),         ### vmin = drift rate in first condition \in (0,\infty)]
    #                          vmax = c( 1, 2.4, 3.8, 5, 8, 12),          ### vmax = mean drift rate in last condition \in (\vmin,\infty)]
    #                          mu_d1 = c(-0.5, 0, 0.5),
    #                          mu_d2 = c(-0.5, 0, 0.5),
    #                          s_v1= c(0.5, 1.4, 3),
    #                          s_v2= c(0.5, 1.4, 3),
    #                          s_d1= c(0.5, 1.4, 3),
    #                          s_d2= c(0.5, 1.4, 3),
    #                          rho_d= c(-0.5, 0, 0.5),
    #                          rho_v= c(-0.5, 0, 0.5),
    #                          t0 = c(0.2, 0.7), ### t0 = minimal non-decision time (proportional to minimum observed RT)
    #                          st0 = seq(0.1, 1.2, length.out=3))    ### st0 = range of (uniform dist) non-decision time
  }
  # Remove columns for fixed parameters
  init_grid <- init_grid[setdiff(names(init_grid), names(fixed))]
  init_grid <- unique(init_grid)

  # Remove rows for which maximum discriminability is lower than minimum discriminability
  init_grid <- init_grid[init_grid$vmin<init_grid$vmax, ]


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

  ## Guess suitable confidence thresholds from theoretical distribution of
  ## the confidence measure and proportion of ratings in the data
  init_thetas <- get_thetas_for_init_grid_MTLNR_simulations(init_grid, df, nRatings, fixed)





  #### 1.1. For Nelder-Mead transform all parameters to real values ####
  if (optim_method=="Nelder-Mead") {
    ## change parametrisation (should be on the whole real line) and
    #  span V-parameters between vmin and vmax equidistantly for all conditions
    inits <- data.frame(matrix(data=NA, nrow= nrow(init_grid),
                               ncol = nConds))
    for (i in 1:nConds){
      inits[,i] <- init_grid[[paste("v", i, sep="")]]
    }

    if (!("mu_d1" %in% names(fixed))) inits <- cbind(inits, init_grid$mu_d1)
    if (!("mu_d2" %in% names(fixed))) inits <- cbind(inits, init_grid$mu_d2)
    if (!("s_v1"  %in% names(fixed))) inits <- cbind(inits, log(init_grid$s_v1 ))
    if (!("s_v2"  %in% names(fixed))) inits <- cbind(inits, log(init_grid$s_v2 ))
    if (!("s_d1"  %in% names(fixed))) inits <- cbind(inits, log(init_grid$s_d1 ))
    if (!("s_d2"  %in% names(fixed))) inits <- cbind(inits, log(init_grid$s_d2 ))
    if (!("rho_d" %in% names(fixed))) inits <- cbind(inits, qnorm(init_grid$rho_d))
    if (!("rho_v" %in% names(fixed))) inits <- cbind(inits, qnorm(init_grid$rho_v))
    if (!("t0"    %in% names(fixed))) inits <- cbind(inits, qnorm(init_grid$t0   ))
    if (!("st0"   %in% names(fixed))) inits <- cbind(inits, log(init_grid$st0  ))

    ## Lowest theta has to be greater than 1, all others have to be ascending
    inits <- cbind(inits, log(init_thetas[,1]))
    if (nRatings > 2) {
      inits <- cbind(inits, log(init_thetas[,-1]))
    }
    if (!sym_thetas) {
      inits <- cbind(inits, log(init_thetas[,1]))
      if (nRatings > 2) {
        inits <- cbind(inits, log(init_thetas[,-1]))
      }
      cols_theta <- c('thetaLower1', rep(paste("dthetaLower", 2:(nRatings-1), sep=""), times=nRatings>2),
                      'thetaUpper1', rep(paste("dthetaUpper", 2:(nRatings-1), sep=""), times=nRatings>2))
    } else {
      cols_theta <- c("theta1", rep(paste("dtheta", 2:(nRatings-1), sep=""), times=nRatings>2))
    }

    ##replace all +-Inf with big/tiny numbers
    inits[inits==Inf]<- 1e6
    inits[inits==-Inf]<- -1e6
    parnames <- c(paste("v", 1:nConds, sep=""), 'mu_d1','mu_d2','s_v1','s_v2',
                  's_d1','s_d2' ,'rho_d','rho_v','t0','st0', cols_theta)
    names(inits) <- setdiff(parnames, names(fixed))
  } else {
    ##### 1.2. For box-constraint optimisation algorithm span drifts and thresholds equidistantly
    if (sym_thetas) {
      init_grid["theta1"] <- init_thetas[,1]
      if (nRatings > 2) {
        for (i in 2:(nRatings-1)) {
          init_grid[paste("dtheta", i, sep="")] <- init_thetas[,i]
        }
        cols_theta <- c("theta1", paste("dtheta", 2:(nRatings-1), sep=""))
      } else {
        cols_theta <- c("theta1")
      }
    } else {
      init_grid[c("thetaUpper1", "thetaLower1")] <- init_thetas[,1]
      if (nRatings > 2) {
        for (i in 2:(nRatings-1)) {
          init_grid[paste(c("dthetaUpper", "dthetaLower"), i, sep="")] <-  init_thetas[,i]
        }
        cols_theta <- c('thetaLower1', paste("dthetaLower", 2:(nRatings-1), sep=""),
                        'thetaUpper1', paste("dthetaUpper", 2:(nRatings-1), sep=""))
      } else {
        cols_theta <- c("thetaLower1", "thetaUpper1")
      }
    }
    #parnames <- c('a', 'b', 't0', 'st0', fitted_weights, paste("v", 1:nConds, sep=""), cols_theta)
    parnames <- c(paste("v", 1:nConds, sep=""), 'mu_d1','mu_d2','s_v1','s_v2',
                  's_d1','s_d2' ,'rho_d','rho_v','t0','st0', cols_theta)
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
    clusterExport(cl, c("df","mint0",
                        "nConds","nRatings", "sym_thetas", "fixed"), envir = environment())
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
                   function(p) try(neglikelihood_MTLNR_free(p, df,  nConds, nRatings, fixed, mint0, sym_thetas, precision),
                                   silent=TRUE))
        #stopCluster(cl)
      } else {
        logL <-
          apply(inits, MARGIN = 1,
                function(p) try(neglikelihood_MTLNR_free(p, df, nConds, nRatings, fixed, mint0,  sym_thetas, precision),
                                silent = TRUE))
      }
    } else {
      if (useparallel) {
        logL <-
          parApply(cl, inits, MARGIN=1,
                   function(p) try(neglikelihood_MTLNR_bounded(p, df, nConds, nRatings, fixed, mint0, sym_thetas, precision),
                                   silent=TRUE))
        #stopCluster(cl)
      } else {
        logL <-
          apply(inits, MARGIN = 1,
                function(p) try(neglikelihood_MTLNR_bounded(p, df, nConds, nRatings, fixed, mint0, sym_thetas, precision),
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
    #                   "v(1:nConds)",    mu_d1,mu_d2,s_v1,s_v2,s_d1,s_d2,rho_d,rho_v,t0,st0, thetaLower1, dthetaLower2.., thetaUpper1... (or theta1,...)
    lower_optbound <- c(rep(-Inf, nConds), -Inf,-Inf,    0,   0,   0,  0,  -1,   -1,  0, 0,   rep(rep(0, nRatings-1), 2-as.numeric(sym_thetas)))[!(parnames %in% names(fixed))]
    upper_optbound <- c(rep( Inf, nConds),  Inf, Inf,  Inf, Inf, Inf, Inf,  1,    1,  1, Inf, rep(Inf, (2-as.numeric(sym_thetas))*(nRatings-1)))[!(parnames %in% names(fixed))]
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
        start <- start + rnorm(length(start), sd=pmax(0.001, abs(t(t(start))/20)))
        if (optim_method == "Nelder-Mead") {
          try(m <- optim(par = start,
                         fn = neglikelihood_MTLNR_free,
                         data=df, nConds=nConds, nRatings=nRatings,
                         fixed=fixed, mint0=mint0,
                         sym_thetas=sym_thetas, precision=precision,
                         method="Nelder-Mead",
                         control = list(maxit = opts$maxit, reltol = opts$reltol)))
        } else if (optim_method =="bobyqa") {
          start <- pmax(pmin(start, upper_optbound-1e-6), lower_optbound+1e-6)
          try(m <- bobyqa(par = start,
                          fn = neglikelihood_MTLNR_bounded,
                          lower = lower_optbound, upper = upper_optbound,
                          data=df, nConds=nConds, nRatings=nRatings,
                          fixed=fixed, mint0=mint0,
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
                         fn = neglikelihood_MTLNR_bounded,
                         lower = lower_optbound, upper = upper_optbound,
                         data=df, nConds=nConds, nRatings=nRatings,
                         fixed=fixed, mint0=mint0,
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
        start <- start + rnorm(length(start), sd=pmax(0.001, abs(t(t(start))/20)))
        if (optim_method == "Nelder-Mead") {
          try(m <- optim(par = start,
                         fn = neglikelihood_MTLNR_free,
                         data=df, nConds=nConds, nRatings=nRatings,
                         fixed=fixed, mint0=mint0,
                         sym_thetas=sym_thetas, precision=precision,
                         method="Nelder-Mead",
                         control = list(maxit = opts$maxit, reltol = opts$reltol)))
        } else if (optim_method =="bobyqa") {
          start <- pmax(pmin(start, upper_optbound-1e-6), lower_optbound+1e-6)
          try(m <- bobyqa(par = start,
                          fn = neglikelihood_MTLNR_bounded,
                          lower = lower_optbound, upper = upper_optbound,
                          data=df, nConds=nConds, nRatings=nRatings,
                          fixed=fixed, mint0=mint0,
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
                         fn = neglikelihood_MTLNR_bounded,
                         lower = lower_optbound, upper = upper_optbound,
                         data=df, nConds=nConds, nRatings=nRatings,
                         fixed=fixed, mint0=mint0,
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
      res[,paste("v",1:(nConds), sep="")] <- p[1:(nConds)]

      if (!("mu_d1" %in% names(fixed))) res$mu_d1 <- p[["mu_d1"]]
      if (!("mu_d2" %in% names(fixed))) res$mu_d2 <- p[["mu_d2"]]
      if (!("s_v1"  %in% names(fixed))) res$s_v1  <- exp(p[["s_v1" ]])
      if (!("s_v2"  %in% names(fixed))) res$s_v2  <- exp(p[["s_v2" ]])
      if (!("s_d1"  %in% names(fixed))) res$s_d1  <- exp(p[["s_d1" ]])
      if (!("s_d2"  %in% names(fixed))) res$s_d2  <- exp(p[["s_d2" ]])
      if (!("rho_d" %in% names(fixed))) res$rho_d <- pnorm(p[["rho_d"]])
      if (!("rho_v" %in% names(fixed))) res$rho_v <- pnorm(p[["rho_v"]])
      if (!("t0"    %in% names(fixed))) res$t0    <- pnorm(p[["t0"   ]])*mint0
      if (!("st0"   %in% names(fixed))) res$st0   <- exp(p[["st0"  ]])

      if (nRatings > 2) {
        if (sym_thetas) {
          res[,paste("theta",1:(nRatings-1), sep="")] <- cumsum(c(exp(p[["theta1"]]), exp(p[paste0("dtheta", 2:(nRatings-1))])))
        } else {
          res[,paste("thetaUpper",1:(nRatings-1), sep="")] <- cumsum(c(exp(p[["thetaUpper1"]]), exp(p[paste0("dthetaUpper", 2:(nRatings-1))])))
          res[,paste("thetaLower",1:(nRatings-1), sep="")] <- cumsum(c(exp(p[["thetaLower1"]]), exp(p[paste0("dthetaLower", 2:(nRatings-1))])))
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
      names(res) <- names(inits)      # a, b, t0, st0, wrt and wint (maybe), v1, v2,....,, thetaLower1,dthetaLower2-4,   thetaUpper1,dthetaUpper2-4,
      if (!("t0" %in% names(fixed))) res$t0 <- res$t0*mint0
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
      # rest by the obvious best-fitting thresholds (e.g. 1/Inf for the lowest/highest...)
      res <- fill_thresholds(res, used_cats, actual_nRatings, 0)
      k <- k+(as.numeric(!sym_thetas)+1)*(actual_nRatings-nRatings)
      nRatings <- actual_nRatings
    }

    if (length(fixed)>=1) {
      res <- cbind(res, as.data.frame(fixed))
    }
    # if (res[['s_d']] == "b") res$a <- res$b
    # if (res[['b']] == "a") res$b <- res$a
    if (sym_thetas) {

      res <- res[,c(paste("v", 1:nConds, sep=""), 'mu_d1','mu_d2','s_v1','s_v2', 's_d1','s_d2' ,'rho_d','rho_v','t0','st0', paste("theta", 1:(nRatings-1), sep=""))]
    } else {
      res <- res[,c(paste("v", 1:nConds, sep=""), 'mu_d1','mu_d2','s_v1','s_v2', 's_d1','s_d2' ,'rho_d','rho_v','t0','st0', paste("thetaLower", 1:(nRatings-1), sep=""), paste("thetaUpper", 1:(nRatings-1), sep=""))]
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


neglikelihood_MTLNR_free <-   function(p, data,
                                     nConds, nRatings,
                                     fixed, mint0, sym_thetas, precision)
{
  # get parameter vector back from real transformations
  paramDf <-  data.frame(matrix(nrow=1, ncol=0))
  if (length(fixed)>=1) {
    paramDf <- cbind(paramDf, as.data.frame(fixed))
  }
  paramDf[,paste("v",1:(nConds), sep="")] <- p[1:(nConds)]
  if (!("mu_d1" %in% names(fixed))) paramDf$mu_d1 <- p[["mu_d1"]]
  if (!("mu_d2" %in% names(fixed))) paramDf$mu_d2 <- p[["mu_d2"]]
  if (!("s_v1"  %in% names(fixed))) paramDf$s_v1  <- exp(p[["s_v1" ]])
  if (!("s_v2"  %in% names(fixed))) paramDf$s_v2  <- exp(p[["s_v2" ]])
  if (!("s_d1"  %in% names(fixed))) paramDf$s_d1  <- exp(p[["s_d1" ]])
  if (!("s_d2"  %in% names(fixed))) paramDf$s_d2  <- exp(p[["s_d2" ]])
  if (!("rho_d" %in% names(fixed))) paramDf$rho_d <- pnorm(p[["rho_d"]])
  if (!("rho_v" %in% names(fixed))) paramDf$rho_v <- pnorm(p[["rho_v"]])
  if (!("t0"    %in% names(fixed))) paramDf$t0    <- pnorm(p[["t0"   ]])*mint0
  if (!("st0"   %in% names(fixed))) paramDf$st0   <- exp(p[["st0"  ]])
  # if (paramDf$a == "b") paramDf$a <- paramDf$b
  # if (paramDf$b == "a") paramDf$b <- paramDf$a

  if (nRatings > 2) {
    if (sym_thetas) {
      paramDf[,paste("theta",1:(nRatings-1), sep="")] <- cumsum(c(exp(p[["theta1"]]), exp(p[paste0("dtheta", 2:(nRatings-1))])))
    } else {
      paramDf[,paste("thetaUpper",1:(nRatings-1), sep="")] <- cumsum(c(exp(p[["thetaUpper1"]]), exp(p[paste0("dthetaUpper", 2:(nRatings-1))])))
      paramDf[,paste("thetaLower",1:(nRatings-1), sep="")] <- cumsum(c(exp(p[["thetaLower1"]]), exp(p[paste0("dthetaLower", 2:(nRatings-1))])))
    }
  } else {
    if (sym_thetas) {
      paramDf[,paste("theta",1:(nRatings-1), sep="")] <- exp(p[["theta1"]])
    } else {
      paramDf[,paste("thetaUpper",1:(nRatings-1), sep="")] <- exp(p[["thetaUpper1"]])
      paramDf[,paste("thetaLower",1:(nRatings-1), sep="")] <- exp(p[["thetaLower1"]])
    }
  }
  if (any(is.infinite(t(paramDf))) || any(is.na(t(paramDf)))) {
    return(1e12)
  }
  negloglik <- -LogLikMTLNR(data, paramDf, precision)

  return(negloglik)
}




neglikelihood_MTLNR_bounded <-   function(p, data,
                                        nConds, nRatings,
                                        fixed, mint0, sym_thetas, precision)
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
      # tryCatch( paramDf[paste("thetaUpper", 2:(nRatings-1), sep="")] <- c(t(paramDf['thetaUpper1'])) + cumsum(c(t(paramDf[grep(names(paramDf), pattern = "dthetaUpper", value=TRUE)]))),
      # error = function(e) {
      #   print(dput(paramDf))
      #   print(names(paramDf))
      #   print(names(p))
      #   print(p)
      #   })
      paramDf[paste("thetaLower", 2:(nRatings-1), sep="")] <- c(t(paramDf['thetaLower1'])) + cumsum(c(t(paramDf[grep(names(paramDf), pattern = "dthetaLower", value=TRUE)])))
    }
    paramDf <- paramDf[ -grep(names(paramDf), pattern="dtheta")]
  }
  if (length(fixed)>=1) {
    paramDf <- cbind(paramDf, as.data.frame(fixed))
  }
  # if (paramDf$a == "b") paramDf$a <- paramDf$b
  # if (paramDf$b == "a") paramDf$b <- paramDf$a
  if (!("t0" %in% names(fixed))) paramDf['t0'] <- paramDf['t0']*mint0
  negloglik <- -LogLikMTLNR(data, paramDf, precision)
  return(negloglik)
}

