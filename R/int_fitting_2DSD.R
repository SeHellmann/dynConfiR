fitting2DSD <- function(df, nConds, nRatings, fixed, sym_thetas,
                        grid_search, init_grid=NULL, optim_method, opts,
                        logging, filename,
                        useparallel, n.cores,
                        restr_tau, precision,
                        used_cats, actual_nRatings){
  ## Be sure that the parallel cluster is stopped if anything happens (error or user interrupt)
  on.exit(try(stopCluster(cl), silent = TRUE))
  ## Set restrictions on tau
  mint0 <-  min(df$rt)
  #if (mint0 < 0.2) {mint0 <- 0.2}
  simult_conf <- FALSE
  if (restr_tau == "simult_conf") {
    # If choice and confidence judgments are given simultaneously, we
    # assume that all judgment processes run sequentially and
    # response time = decision time + interjudgment time (tau)
    #                  + non-judgment time (t0)
    simult_conf = TRUE
    restr_tau = 1
  }
  ### 1. Generate initial grid for grid search over possible parameter sets ####
  #### Create grid ####
  if (is.null(init_grid)) {
    if (restr_tau == Inf) {
      tau = seq(0.3, 1.5, length.out = 3)
    } else if (simult_conf) {
      tau = seq(0.1, 0.8, length.out = 3)
    } else {
      if (!(is.numeric(restr_tau) && restr_tau >0)) {stop(paste("restr_tau must be numeric and positive, Inf or 'simult_conf'. But restr_tau=", restr_tau, sep=""))}
      tau = seq(0.3*restr_tau,0.9*restr_tau, length.out = 3)
    }
    init_grid <- expand.grid(a = c(1.3,2, 2.8, 5),
                             vmin = c(0.01, 0.1, 1.3),
                             vmax = c(1.4, 2.5, 3.7, 5),
                             sv = c(0.01, 0.8, 1.5),
                             z = sum((df$response==1)*df$n)/sum(df$n),
                             sz = c(0.01, 0.4),
                             t0 = c(0.4, 0.8),
                             st0 = c(0.1, 0.3, 0.7),
                             tau = tau,
                             lambda=c(0, 0.5, 1, 2))
  }

  # Remove columns for fixed parameters
  init_grid <- init_grid[setdiff(names(init_grid), names(fixed))]
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
  ## Guess suitable confidence thresholds from theoretical distribution of
  ## the confidence measure and proportion of ratings in the data
  init_thetas <- get_thetas_for_init_grid_dynWEV_simulations(init_grid, df, nRatings, simult_conf,
                                                 c(fixed, w=1, sigvis=1, svis=1), # 2DSD is dynWEV with w fixed to 1
                                                 mint0 = mint0)


  #### For Nelder-Mead transform all parameters to real values ####
  if (optim_method=="Nelder-Mead") {
    ## change parametrization (should be on the whole real line) and
    inits <- data.frame(matrix(data=NA, nrow= nrow(init_grid),
                               ncol = nConds))
    for (i in 1:nConds){
        inits[,i] <- log(init_grid[[paste("v", i, sep="")]])
    }

    if (!("a" %in% names(fixed))) inits <- cbind(inits, log(init_grid$a))
    if (!("sv" %in% names(fixed))) inits <-  cbind(inits, log(init_grid$sv)) # one SV (SD of drift rate) for all the different conditions
    if (!("z" %in% names(fixed))) inits <- cbind(inits, qnorm(init_grid$z))
    if (!("sz" %in% names(fixed))) inits <- cbind(inits, qnorm(init_grid$sz))
    if (!("t0" %in% names(fixed))) inits <- cbind(inits, qnorm(init_grid$t0))
    if (!("st0" %in% names(fixed))) inits <- cbind(inits, log(init_grid$st0))
    if (!("tau" %in% names(fixed))) {
      if (restr_tau == Inf) {
        inits <- cbind(inits, log(init_grid$tau))
      } else {
        inits <- cbind(inits, qnorm(init_grid$tau / restr_tau))
      }
    }
    if (!("lambda" %in% names(fixed))) inits <- cbind(inits, log(init_grid$lambda))
    inits <- cbind(inits, init_thetas[,1])
    if (nRatings > 2) {
      inits <- cbind(inits, log(init_thetas[,-1]))
    }
    if (!sym_thetas) {
      inits <- cbind(inits, init_thetas[,1])
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
    parnames <- c(paste("v", 1:nConds, sep=""), 'a', 'sv', 'z', 'sz', 't0', 'st0', 'tau', 'lambda', cols_theta)
    names(inits) <- setdiff(parnames, names(fixed))

  } else {
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
    parnames <- c('a', 'z', 'sz', paste("v", 1:nConds, sep=""),
                   'st0', 'sv', 't0', cols_theta,'tau', 'lambda')
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
    clusterExport(cl, c("df", "restr_tau", "fixed", "mint0",
                        "nConds","nRatings", "sym_thetas", "simult_conf", "precision"), envir = environment())
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
                   function(p) try(neglikelihood_2DSD_free(p, df,
                                                           restr_tau, nConds, nRatings, fixed, mint0, simult_conf, sym_thetas, precision),
                                   silent=TRUE))
        #stopCluster(cl)
      } else {
        logL <-
          apply(inits, MARGIN = 1,
                function(p) try(neglikelihood_2DSD_free(p, df, restr_tau, nConds, nRatings, fixed, mint0, simult_conf, sym_thetas, precision),
                                silent = TRUE))
      }
    } else {
      if (useparallel) {
        logL <-
          parApply(cl, inits, MARGIN=1,
                   function(p) try(neglikelihood_2DSD_bounded(p, df, restr_tau, nConds, nRatings, fixed, mint0, simult_conf, sym_thetas, precision),
                                   silent=TRUE))
        #stopCluster(cl)
      } else {
        logL <-
          apply(inits, MARGIN = 1,
                function(p) try(neglikelihood_2DSD_bounded(p, df, restr_tau, nConds, nRatings, fixed, mint0, simult_conf,sym_thetas, precision),
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
                      # a,  z, sz,  v1, v2,....,,   st0, sv, t0,thetaLower1, dthetaLower2.., thetaUpper1... (or theta1,...),    tau,      lambda
    lower_optbound <- c(0,  0,  0,  rep(0, nConds), 0,   0,  0, rep(c(-Inf,     rep(0, nRatings-2)), 2-as.numeric(sym_thetas)), 0,        0)[!(parnames %in% names(fixed))]
    upper_optbound <- c(Inf,1,  1,  rep(Inf,nConds),Inf, Inf,1, rep(Inf, (2-as.numeric(sym_thetas))*(nRatings-1)),              restr_tau,Inf)[!(parnames %in% names(fixed))]
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
          m <- try(optim(par = start,
                         fn = neglikelihood_2DSD_free,
                         data=df, restr_tau = restr_tau, nConds=nConds, nRatings=nRatings,
                         fixed = fixed,  mint0 =mint0,
                         simult_conf = simult_conf,sym_thetas=sym_thetas, precision=precision,
                         method="Nelder-Mead",
                         control = list(maxit = opts$maxit, reltol = opts$reltol)))
        } else if (optim_method =="bobyqa") {
          start <- pmax(pmin(start, upper_optbound-1e-6), lower_optbound+1e-6)
          m <- try(bobyqa(par = start,
                          fn = neglikelihood_2DSD_bounded,
                          lower = lower_optbound, upper = upper_optbound,
                          data=df, restr_tau = restr_tau, nConds=nConds, nRatings=nRatings,
                          fixed = fixed,  mint0 =mint0,
                          simult_conf = simult_conf,sym_thetas=sym_thetas, precision=precision,
                          control = list(maxfun=opts$maxfun,
                                         rhobeg = min(0.2, 0.2*restr_tau, 0.2*max(abs(start))),
                                         npt = length(start)+5)))
          ## rhobeg should be: about 0.1*(greatest expected change in parameters --> <= 1-2 (for a, thetas or v's) )
          ##                   smaller than min(abs(upper-lower)) = min(1, restr_tau)
          ##                   --> so we use min(0.2*restr_tau, 0.2, 0.2*max(abs(par))).
          ##                   Default would be: min(0.95, 0.2*max(abs(par))), respectively 0.2*max(upper_optbound-lower_optbound)
          ## rhoend: use default of 1e-6*rhobeg
          if (exists("m") && !inherits(m, "try-error") && !is.na(m$fval)){

            m$value <- m$fval
          }
        } else if (optim_method=="L-BFGS-B") {  ### ToDo: use dfoptim or pracma::grad as gradient!
          start <- pmax(pmin(start, upper_optbound-1e-6), lower_optbound+1e-6)
          m <- try(optim(par = start,
                         fn = neglikelihood_2DSD_bounded,
                         lower = lower_optbound, upper = upper_optbound,
                         data=df,  restr_tau = restr_tau, nConds=nConds, nRatings=nRatings,
                         fixed = fixed,  mint0 =mint0,
                         simult_conf = simult_conf,sym_thetas=sym_thetas, precision=precision,
                         method="L-BFGS-B",
                         control = list(maxit = opts$maxit, factr = opts$factr)))
        } else {
          stop(paste("Not implemented or unknown method: ", optim_method, ". Use 'bobyqa', Nelder-Mead' or 'L-BFGS-B' instead.", sep=""))
        }
        if (logging==TRUE) {
          logger::log_info(paste("Finished attempt No.", i, " restart no. ", l))
        }
        if (!exists("m") || inherits(m, "try-error") || is.na(m$fval)){
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
      for (l in 1:opts$nRestarts){
        start <- start + rnorm(length(start), sd=pmax(0.001, abs(t(t(start))/20)))
        names(start) <- parnames
        if (optim_method == "Nelder-Mead") {
          m <- try(optim(par = start,
                         fn = neglikelihood_2DSD_free,
                         data=df,  restr_tau = restr_tau, nConds=nConds, nRatings=nRatings,
                         fixed = fixed,  mint0 =mint0,
                         simult_conf = simult_conf,sym_thetas=sym_thetas, precision=precision,
                         method="Nelder-Mead",
                         control = list(maxit = opts$maxit, reltol = opts$reltol)))
        } else if (optim_method =="bobyqa") {
          start <- pmax(pmin(start, upper_optbound-1e-6), lower_optbound+1e-6)
          m <- try(bobyqa(par = start,
                          fn = neglikelihood_2DSD_bounded,
                          lower = lower_optbound, upper = upper_optbound,
                          data=df,  restr_tau = restr_tau, nConds=nConds, nRatings=nRatings,
                          fixed = fixed, mint0 =mint0,
                          simult_conf = simult_conf,sym_thetas=sym_thetas, precision=precision,
                          control = list(maxfun=opts$maxfun,
                                         rhobeg = min(0.2, 0.2*restr_tau, 0.2*max(abs(start))),
                                         npt = length(start)+5)),
              silent=TRUE)
          ## rhobeg should be: about 0.1*(greatest expected change in parameters --> <= 1-2 (for a, thetas or v's) )
          ##                   smaller than min(abs(upper-lower)) = min(1, restr_tau)
          ##                   --> so we use min(0.2*restr_tau, 0.2, 0.2*max(abs(par))).
          ##                   Default would be: min(0.95, 0.2*max(abs(par))), respectively 0.2*max(upper_optbound-lower_optbound)
          ## rhoend: use default of 1e-6*rhobeg
          if (exists("m") && !inherits(m, "try-error") && !is.na(m$fval)){
            m$value <- m$fval
          }
        } else if (optim_method=="L-BFGS-B") {  ### ToDo: use dfoptim or pracma::grad as gradient!
          start <- pmax(pmin(start, upper_optbound-1e-6), lower_optbound+1e-6)
          m <- try(optim(par = start,
                         fn = neglikelihood_2DSD_bounded,
                         lower = lower_optbound, upper = upper_optbound,
                         data=df,  restr_tau = restr_tau, nConds=nConds, nRatings=nRatings,
                         fixed = fixed,  mint0 =mint0,
                         simult_conf = simult_conf,sym_thetas=sym_thetas, precision=precision,
                         method="L-BFGS-B",
                         control = list(maxit = opts$maxit, factr = opts$factr)))
        } else {
          stop(paste("Not implemented or unknown method: ", optim_method, ". Use 'bobyqa', Nelder-Mead' or 'L-BFGS-B' instead.", sep=""))
        }
        if (!exists("m") || inherits(m, "try-error") || is.na(m$fval)){
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
        return(rep(NA, length(start)+1))
      } # end of node-function
    }
    clusterExport(cl, c("parnames", "opts", "optim_method","optim_node"), envir = environment())
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
    # Save the number of actually fitted parameters by the optimization
    k <- length(fit$par)
    N <- nrow(df)
    p <- fit$par

    if (optim_method=="Nelder-Mead") {
      res[,paste("v",1:(nConds), sep="")] <- exp(p[1:(nConds)])
      if (length(fixed)>=1) res <- cbind(res, as.data.frame(fixed))

      if (!("a" %in% names(fixed))) res$a <- exp(p[["a"]])
      if (!("sv" %in% names(fixed))) res$sv <- exp(p[["sv"]])
      if (!("z" %in% names(fixed))) res$z <- pnorm(p[["z"]]) ## relative mean starting point
      if (!("sz" %in% names(fixed))) res$sz <- (min(res$z, (1-res$z))*2)*pnorm(p[["sz"]])
      if (!("t0" %in% names(fixed))) res$t0 <- pnorm(p[["t0"]])*mint0
      if (!("st0" %in% names(fixed))) res$st0 <- exp(p[["st0"]])
      if (!("tau" %in% names(fixed))) {
        if (restr_tau == Inf) {
          res$tau <- exp(p[["tau"]])
        } else if (simult_conf) {
          res$tau <- pnorm(p[["tau"]])*(mint0-res$t0)
        } else {
          res$tau <- restr_tau * pnorm(p[["tau"]])
        }
      }
      if (!("lambda" %in% names(fixed))) res$lambda <- exp(p[["lambda"]])
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
      p <- c(t(p))
      names(p) <- names(inits)
      if (nRatings>2) {
        if (sym_thetas) {
          p[paste("theta", 2:(nRatings-1), sep="")] <- c(t(p['theta1'])) + c(t(cumsum(p[grep(names(p), pattern = "dtheta", value=TRUE)])))
        } else {
          p[paste("thetaUpper", 2:(nRatings-1), sep="")] <- c(t(p['thetaUpper1'])) + cumsum(c(t(p[grep(names(p), pattern = "dthetaUpper", value=TRUE)])))
          p[paste("thetaLower", 2:(nRatings-1), sep="")] <- c(t(p['thetaLower1'])) + cumsum(c(t(p[grep(names(p), pattern = "dthetaLower", value=TRUE)])))
        }
        p <- p[ -grep(names(p), pattern="dtheta")]
      }
      if (length(fixed)>=1) p <- c(p, unlist(fixed))
      if (!("t0" %in% names(fixed))) p['t0'] <- p['t0']*mint0
      if (simult_conf & !("tau" %in% names(fixed))) p['tau'] <- p['tau']*(mint0-p[['t0']])
      if (!("sz" %in% names(fixed))) p['sz'] <- (min(p['z'], (1-p['z']))*2)*p["sz"] ##

      res <-   data.frame(matrix(nrow=1, ncol=length(p)))
      res[1,] <- p
      names(res) <- names(p)      # a,  z, sz,v1, v2,....,,   st0, sv, t0, thetaLower1,dthetaLower2-4,   thetaUpper1,dthetaUpper2-4,    tau

    }
    if (!is.null(used_cats)) {
      # If some rating categories are not used, we fit less thresholds numerically and fill up the
      # rest by the obvious best-fitting thresholds (e.g. +/- Inf for the lowest/highest...)
      res <- fill_thresholds(res, used_cats, actual_nRatings, -1e+24)
      # Get number of truely fitted parameters by adding the added confidence thresholds
      k <- k + (as.numeric(!sym_thetas)+1)*(actual_nRatings-nRatings)
      nRatings <- actual_nRatings
    }

    if (sym_thetas) {
      parnames <- c(paste("v", 1:nConds, sep=""), 'sv', 'a', 'z', 'sz', 't0','st0', paste("theta", 1:(nRatings-1), sep=""), 'tau', 'lambda')
    } else {
      parnames <- c(paste("v", 1:nConds, sep=""), 'sv', 'a', 'z', 'sz', 't0','st0', paste("thetaLower", 1:(nRatings-1), sep=""), paste("thetaUpper", 1:(nRatings-1), sep=""), 'tau', 'lambda')
    }
    res <- res[, parnames]

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




neglikelihood_2DSD_free <-   function(p, data,
                                      restr_tau, nConds, nRatings, fixed, mint0, simult_conf, sym_thetas, precision=3)
{
  # get parameter vector back from real transformations
  paramDf <-  data.frame(matrix(nrow=1, ncol=0))
  paramDf[,paste("v",1:(nConds), sep="")] <- exp(p[1:(nConds)])
  if (length(fixed)>=1) paramDf <- cbind(paramDf, as.data.frame(fixed))

  if (!("a" %in% names(fixed))) paramDf$a <- exp(p[["a"]])
  if (!("sv" %in% names(fixed))) paramDf$sv <- exp(p[["sv"]])
  if (!("z" %in% names(fixed))) paramDf$z <- pnorm(p[["z"]]) ## relative mean starting point
  if (!("sz" %in% names(fixed))) paramDf$sz <- (min(paramDf$z, (1-paramDf$z))*2)*pnorm(p[["sz"]]) ##
  if (!("t0" %in% names(fixed))) paramDf$t0 <- pnorm(p[["t0"]])*mint0
  if (!("st0" %in% names(fixed))) paramDf$st0 <- exp(p[["st0"]])
  if (!("tau" %in% names(fixed))) {
    if (restr_tau == Inf) {
      paramDf$tau <- exp(p[["tau"]])
    } else if (simult_conf) {
      paramDf$tau <- pnorm(p[["tau"]])*(mint0-paramDf$t0)
    } else {
      paramDf$tau <- restr_tau * pnorm(p[["tau"]])
    }
  }
  if (!("lambda" %in% names(fixed))) paramDf$lambda <- exp(p[["lambda"]])
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
  negloglik <- -LogLikWEV(data, paramDf, "2DSD", simult_conf, precision, stop_on_error=FALSE)
  return(negloglik)
}


neglikelihood_2DSD_bounded <-   function(p, data,
                                         restr_tau, nConds, nRatings, fixed, simult_conf, mint0, sym_thetas=FALSE, precision=3)
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
  if (!("t0" %in% names(fixed))) paramDf['t0'] <- paramDf['t0']*mint0
  if (!("sz" %in% names(fixed))) paramDf['sz'] <- (min(paramDf['z'], 1-paramDf['z'])*2)*paramDf['sz']
  if (simult_conf & !("tau" %in% names(fixed))) {
    paramDf['tau'] <- paramDf['tau']*(mint0-paramDf['t0'])
  }
  #paramDf[grep(names(paramDf), pattern = "thetaLower", value=TRUE)] <- paramDf$a - paramDf[grep(names(paramDf), pattern = "thetaLower", value=TRUE)]
  negloglik <- -LogLikWEV(data, paramDf, "2DSD", simult_conf, precision, stop_on_error=FALSE)
  return(negloglik)
}
