# fittingdynWEV_formula <- function(DVs, model_matrix,
#                                   fixed,
#                                   nRatings, sym_thetas,
#                           optim_method, opts,
#                           logging, filename,
#                           useparallel, n.cores,
#                           restr_tau, precision,
#                           used_cats, actual_nRatings){
#
#   model_matrix
#   fixed
#   fit_pars_beta_names
#   const_pars
#
#   #fit_par_names
#   fit_pars_columns
#
#
#   #manipulations
#
#   #fit_pars_betanames <- sub("\\(Intercept\\)", "1",fit_pars_columns_names)
#
#   par_thetas <- if (!sym_thetas) {
#     if (nRatings>2) {
#       c("thetaLower1", "thetaUpper1", paste0("dtheta", c("Lower", "Upper"), rep(2:(nRatings-1), each=2)))
#     } else {
#       c("thetaLower", "thetaUpper")
#     }
#   } else {
#     if (nRatings>2) {
#       c("theta1", paste0("dtheta", rep(2:(nRatings-1), each=2)))
#     } else {
#       c("theta")
#     }
#   }
#
#
#
#
#
#   ## Be sure that the parallel cluster is stopped if anything happens (error or user interupt)
#   on.exit(try(stopCluster(cl), silent = TRUE))
#   ## Set restrictions on tau
#   maxt0 <-  min(DVs$rt)
#   simult_conf <- FALSE
#   if (restr_tau == "simult_conf") {
#     # If choice and confidence judgments are given simultaneously, we
#     # assume that all judgment processes ran sequentially and
#     # response time = decision time + interjudgment time (tau)
#     #                  + non-judgment time (t0)
#     simult_conf = TRUE
#     restr_tau = 1
#   }
#   ### 1. Generate initial grid for grid search over possible parameter sets ####
#   #### Create grid ####
#   if (!simult_conf && !(is.numeric(restr_tau) && restr_tau >0)) {stop(paste("restr_tau must be numeric and positive, Inf or 'simult_conf'. But restr_tau=", restr_tau, sep=""))}
#
#
#   #beta <- rnorm(length(fit_pars_columns_names)+length(const_pars)+length(par_thetas))
#   #names(beta) <-
#   beta_names <- c(fit_pars_beta_names, const_pars, par_thetas)
#   k <- length(beta_names)
#
#   n_initials <- 100
#   inits <- matrix(runif(n_initials*k), nrow = n_initials, dimnames = list(NULL, beta_names))
#   if ("st0" %in% beta_names) inits[,c("st0")] <- inits[,c("st0")]/2 - 1.6  # Rescale st0 initials to lower values, because integration would otherwise take a lot of time
#
#   ### If no grid-search is desired use mean of possible parameters
#   if (!grid_search) {
#     inits <- colMeans(inits)
#   }
#
#
#   ## Intermezzo: Setup cluster for parallelization   ####
#   if (useparallel) {
#     if (is.null(n.cores)) {
#       n.cores <- detectCores()-1
#     }
#     cl <- makeCluster(type="SOCK", n.cores)
#     clusterExport(cl, c("DVs", "restr_tau", "maxt0",
#                         "nConds","nRatings", "fixed", "simult_conf", "sym_thetas", "precision"), envir = environment())
#   }
#
#
#   ### 2. Search initial grid before optimization  ####
#   if (grid_search) {
#     if (logging==TRUE) {
#       logger::log_info(paste(length(inits[,1]), "...parameter sets to check"))
#       logger::log_info(paste("data got ", nrow(DVs), " rows"))
#       t00 <- Sys.time()
#       logger::log_info("Searching initial values ...")
#     }
#     t0 <- Sys.time()
#     if (useparallel) {
#         logL <-
#           parApply(cl, inits, MARGIN=1,
#                    function(p) try(neglikelihood_formula(p, DVs,
#                                                          model_matrix, fit_pars_columns,
#                                                          fixed, restr_tau, nRatings, maxt0, simult_conf, sym_thetas, precision),
#                                    silent=TRUE))
#         #stopCluster(cl)
#       } else {
#         logL <-
#           apply(inits, MARGIN = 1,
#                 function(p) try(neglikelihood_formula(p, DVs,
#                                                       model_matrix, fit_pars_columns,
#                                                       fixed, restr_tau, nRatings, maxt0, simult_conf, sym_thetas, precision),
#                                 silent = TRUE))
#       }
#   print(logL)
#   t2 <- difftime(Sys.time(), t0, units = "min")
#   logL1 <- logL
#     logL <- as.numeric(logL)
#
#     inits <- inits[order(logL),]
#     if (logging==TRUE) {
#       logger::log_success(paste("Initial grid search took...",as.character(round(as.double(difftime(Sys.time(),t00,units = "mins")), 2))," mins"))
#     }
#   } else {
#     logL <- NULL
#   }
#
#   #### 3. Optimization ####
#   if (logging==TRUE) {
#     logger::log_info("Start fitting ... ")
#   }
#   if (!useparallel || (opts$nAttempts==1)) {
#     noFitYet <- TRUE
#     for (i in 1:opts$nAttempts){
#       start <- c(t(inits[i,]))
#       names(start) <- colnames(inits)
#       for (l in 1:opts$nRestarts){
#         start <- start + rnorm(length(start), sd=pmax(0.001, abs(t(t(start))/20)))
#         if (optim_method == "Nelder-Mead") {
#           try(m <- optim(par = start,
#                          fn = neglikelihood_formula,
#                          DVs=DVs,
#                          model_matrix = model_matrix,
#                          fit_pars_columns=fit_pars_columns,
#                          restr_tau = restr_tau, nRatings=nRatings,
#                          fixed=fixed, maxt0=maxt0, simult_conf=simult_conf,
#                          sym_thetas=sym_thetas, precision=precision,
#                          method="Nelder-Mead",
#                          control = list(maxit = opts$maxit, reltol = opts$reltol)))
#         } else if (optim_method =="bobyqa") {
#           try(m <- bobyqa(par = start, lower=-Inf, upper=Inf,
#                           fn = neglikelihood_formula,
#                           DVs=DVs,
#                           model_matrix = model_matrix,
#                           fit_pars_columns=fit_pars_columns,
#                           restr_tau = restr_tau, nRatings=nRatings,
#                           fixed=fixed, maxt0=maxt0, simult_conf=simult_conf,
#                           sym_thetas=sym_thetas, precision=precision,
#                           control = list(maxfun=opts$maxfun,
#                                          rhobeg = 2,
#                                          npt = length(start)+5)))
#           ## rhobeg should be: about 0.1*(greatest expected change in parameters --> <= 1-2 (for a, thetas or v's) )
#           ##                   smaller than min(abs(upper-lower)) = min(1, restr_tau)
#           ##                   --> so we use min(0.2*restr_tau, 0.2, 0.2*max(abs(par))).
#           ##                   Default would be: min(0.95, 0.2*max(abs(par))), respectively 0.2*max(upper_optbound-lower_optbound)
#           ## rhoend: use default of 1e-6*rhobeg
#           if (exists("m") && !inherits(m, "try-error")){
#             m$value <- m$fval
#           }
#         } else {
#           stop(paste("Not implemented or unknown method: ", optim_method, ". Use 'bobyqa'or 'Nelder-Mead' instead.", sep=""))
#         }
#         if (logging==TRUE) {
#           logger::log_info(paste("Finished attempt No.", i, " restart no. ", l))
#         }
#         if (!exists("m") || inherits(m, "try-error")){
#           if (logging==TRUE) {
#             logger::log_error(paste("No fit obtained at attempt No.", i))
#             logger::log_error(paste("Used parameter set", paste(start, sep="", collapse=" "), sep=" ", collapse = ""))
#           }
#           break
#         }
#         if (exists("m") && is.list(m)){
#           if (noFitYet) {
#             fit <- m
#             noFitYet <- FALSE
#             if (logging==TRUE) {
#               logger::log_info(paste("First fit obtained at attempt No.", i))
#               attempt <- i
#               save(logL, inits,  DVs,fit, attempt,file=filename)
#             }
#             start <- fit$par
#             names(start) <- names(inits)
#           } else if (m$value < fit$value) {
#             fit <- m
#             if (logging==TRUE) {
#               logger::log_info(paste("New fit at attempt No.", i, " restart no. ", l))
#               attempt <- i
#               save(logL, inits,  DVs,fit, attempt,file=filename)
#             }
#             start <- fit$par
#             names(start) <- names(inits)
#           } # end of if better value
#         }   # end of if we got a optim-result at all
#       }     # end of for restarts
#     }       # end of for initial start values
#   } else {  # if useparallel
#     starts <- inits[(1:opts$nAttempts),]
#     parnames <- names(starts)
#
#     optim_node <- function(start) { # define optim-routine to run on each node
#       noFitYet <- TRUE
#       start <- c(t(start))
#       for (l in 1:opts$nRestarts){
#         start <- start + rnorm(length(start), sd=pmax(0.001, abs(t(t(start))/20)))
#         names(start) <- parnames
#         if (optim_method == "Nelder-Mead") {
#           m <- try(optim(par = start,
#                          fn = neglikelihood_formula,
#                          data=DVs,  restr_tau = restr_tau, nConds=nConds, nRatings=nRatings,
#                          fixed=fixed, maxt0=maxt0, simult_conf=simult_conf,
#                          sym_thetas=sym_thetas, precision=precision,
#                          method="Nelder-Mead",
#                          control = list(maxit = opts$maxit, reltol = opts$reltol)))
#         } else if (optim_method =="bobyqa") {
#           m <- try(bobyqa(par = start,
#                           fn = neglikelihood_dynWEV_bounded,
#                           lower = lower_optbound, upper = upper_optbound,
#                           data=DVs,  restr_tau = restr_tau, nConds=nConds, nRatings=nRatings,
#                           fixed=fixed, maxt0=maxt0, simult_conf=simult_conf,
#                           sym_thetas=sym_thetas, precision=precision,
#                           control = list(maxfun=opts$maxfun,
#                                          rhobeg = min(0.2, 0.2*restr_tau, 0.2*max(abs(start))),
#                                          npt = length(start)+5)))
#           ## rhobeg should be: about 0.1*(greatest expected change in parameters --> <= 1-2 (for a, thetas or v's) )
#           ##                   smaller than min(abs(upper-lower)) = min(1, restr_tau)
#           ##                   --> so we use min(0.2*restr_tau, 0.2, 0.2*max(abs(par))).
#           ##                   Default would be: min(0.95, 0.2*max(abs(par))), respectively 0.2*max(upper_optbound-lower_optbound)
#           ## rhoend: use default of 1e-6*rhobeg
#           if (exists("m") && !inherits(m, "try-error")){
#             m$value <- m$fval
#           }
#         } else {
#           stop(paste("Not implemented or unknown method: ", optim_method, ". Use 'bobyqa' or 'Nelder-Mead' instead.", sep=""))
#         }
#         if (!exists("m") || inherits(m, "try-error")){
#           break
#         }
#         if (exists("m") && is.list(m)){
#           if (noFitYet) {
#             fit <- m
#             noFitYet <- FALSE
#             start <- fit$par
#             names(start) <- parnames
#           } else if (m$value < fit$value) {
#             fit <- m
#             start <- fit$par
#             names(start) <- parnames
#           }
#         }
#       }
#       if (exists("fit") && is.list(fit)){
#         return(c(fit$value,fit$par))
#       } else {
#         return(c(m[1] , rep(NA, length(start))))
#       } # end of node-function
#     }
#     clusterExport(cl, c("parnames", "opts", "optim_method","optim_node" ), envir = environment())
#     if (optim_method!="Nelder-Mead") {
#       clusterExport(cl, c("lower_optbound", "upper_optbound"), envir = environment())
#     }
#     optim_outs <- parApply(cl, starts,MARGIN=1, optim_node )
#     stopCluster(cl)
#     optim_outs <- t(optim_outs)
#     best_res <- optim_outs[order(optim_outs[,1]),][1,]
#     fit <- list(par = best_res[-1], value=best_res[1])
#   }   # end of if-else useparallel
#
#   #### 4. Wrap up results ####
#   res <-  list()
#   if(exists("fit") && is.list(fit)){
#     res$k <- length(fit$par)
#     res$N <- nrow(DVs)
#     res$beta <- fit$par
#     names(res$beta) <- names(inits)
#     res$model_matrix <- model_matrix
#     res$manipulations
#
#     res$fixed <- paste(c("sym_thetas", names(fixed)), c(sym_thetas,fixed), sep="=", collapse = ", ")
#     res$negLogLik <- fit$value
#     res$BIC <-  2 * fit$value + res$k * log(res$N)
#     res$AICc <- 2 * fit$value + res$k * 2 + 2*res$k*(res$k-1)/(res$N-res$k-1)
#     res$AIC <- 2 * fit$value + res$k * 2
#
#
#   #   compute parameter vector for res
#       # if (nRatings > 2) {
#       #   if (sym_thetas) {
#       #     res[,paste("theta",1:(nRatings-1), sep="")] <- cumsum(c(p[["theta1"]], exp(p[paste0("dtheta", 2:(nRatings-1))])))
#       #   } else {
#       #     res[,paste("thetaUpper",1:(nRatings-1), sep="")] <- cumsum(c(p[["thetaUpper1"]], exp(p[paste0("dthetaUpper", 2:(nRatings-1))])))
#       #     res[,paste("thetaLower",1:(nRatings-1), sep="")] <- cumsum(c(p[["thetaLower1"]], exp(p[paste0("dthetaLower", 2:(nRatings-1))])))
#       #   }
#       # } else {
#       #   if (sym_thetas) {
#       #     res[,paste("theta",1:(nRatings-1), sep="")] <- p[["theta1"]]
#       #   } else {
#       #     res[,paste("thetaUpper",1:(nRatings-1), sep="")] <- p[["thetaUpper1"]]
#       #     res[,paste("thetaLower",1:(nRatings-1), sep="")] <- p[["thetaLower1"]]
#       #   }
#       # }
#     # p <- fit$par
#     # if (optim_method=="Nelder-Mead") {
#     #
#     # if (!is.null(used_cats)) {
#     #   # If some rating categories are not used, we fit less thresholds numerically and fill up the
#     #   # rest by the obvious best-fitting thresholds (e.g. +/- Inf for the lowest/highest...)
#     #   res <- fill_thresholds(res, used_cats, actual_nRatings, -1e+24)
#     #   nRatings <- actual_nRatings
#     #   k <- ncol(res)
#     # }
#     # if (sym_thetas) {
#     #   parnames <- c(paste("v", 1:nConds, sep=""), 'sv', 'a', 'z', 'sz', 't0','st0', paste("theta", 1:(nRatings-1), sep=""), 'tau', 'w', 'svis','sigvis', 'lambda')
#     # } else {
#     #   parnames <- c(paste("v", 1:nConds, sep=""), 'sv', 'a', 'z', 'sz', 't0','st0', paste("thetaLower", 1:(nRatings-1), sep=""), paste("thetaUpper", 1:(nRatings-1), sep=""), 'tau', 'w', 'svis','sigvis', 'lambda')
#     # }
#     # res <- res[, parnames]
#
#
#     if (logging==TRUE) {
#       logger::log_success("Done fitting and autosaved results")
#       save(logL, DVs, model_matrix, fit, inits, res, file=filename)
#     }
#   }
#   return(res)
# }
#
#
#
# #beta <- inits[1,]
# neglikelihood_formula <-   function(beta, DVs,
#                                     model_matrix, fit_pars_columns,
#                                     fixed,
#                                     restr_tau, nRatings, maxt0, simult_conf, sym_thetas, precision=1e-5)
# {
#   # get parameter vector back from real transformations
#
#   #parnames <- c("v", "z", "a", "sz", "t0", "st0", "d", "sv", "tau", "w", "svis", "sigvis", "lambda", "s", "th1", "th2")
#   parnames <- c("a", "v", "t0", "d","z",  "sz", "sv", "st0","tau", "th1", "th2", "lambda", "w", "muvis", "sigvis", "svis", "s")
#   parammatrix <-  matrix(NA, nrow=nrow(DVs), ncol=length(parnames))
#   colnames(parammatrix) <- parnames
#   fixed01 <- c("z", "sz", "w", "d", "st0")
#   fixedpos <- c("v", "a", "t0", "sv", "svis", "sigvis", "lambda", "s")
#
#   for (i in 1:length(parnames)) {
#     if (parnames[i] %in% names(fixed)) {
#       parammatrix[,parnames[i]] <- fixed[[parnames[i]]]
#     } else {
#       if (parnames[i] %in% names(beta)) parammatrix[,parnames[i]] <- beta[parnames[i]]
#       if (parnames[i] %in% names(fit_pars_columns)) parammatrix[,parnames[i]] <-
#           model_matrix[,fit_pars_columns[[parnames[i]]], drop=FALSE] %*%
#             beta[grepl(names(beta), pattern=paste0(parnames[i], "_"))]
#       if (parnames[i] %in% fixed01) parammatrix[,parnames[i]] <- pnorm(parammatrix[,parnames[i]])
#       if (parnames[i] %in% fixedpos) parammatrix[,parnames[i]] <- exp(parammatrix[,parnames[i]])
#       if (parnames[i] == "sz") parammatrix[,"sz"] <- (pmin(parammatrix[,"z"], (1-parammatrix[,"z"]))*2)*parammatrix[,"sz"]
#       if (parnames[i] == "t0") parammatrix[,"t0"] <- maxt0*parammatrix[,"t0"]
#       if (parnames[i] == "d") parammatrix[,"d"] <- parammatrix[,"t0"]*parammatrix[,"d"]
#       if (parnames[i] == "tau") {
#         if (restr_tau == Inf) {
#           parammatrix[,parnames[i]] <- exp(parammatrix[,parnames[i]])
#         } else if (simult_conf) {
#           parammatrix[,parnames[i]] <- pnorm(parammatrix[,parnames[i]])*(maxt0-parammatrix[,"t0"])
#         } else {
#           parammatrix[,parnames[i]] <- restr_tau * pnorm(parammatrix[,parnames[i]])
#         }
#       }
#     }
#   }
#   #parammatrix[,"v"] <- parammatrix[,"v"]*DVs$stimulus
#   par_thetas <- grep(names(beta), pattern="theta", value=TRUE)
#   if (sym_thetas) {
#       Thetas <- beta[par_thetas]
#       Thetas <- c(-1e+21, Thetas[1] + c(0, cumsum(exp(Thetas[-1]))), Inf)
#       Thetas <- c(Thetas, Thetas)
#     } else {
#       ThetasLower <- beta[grep(par_thetas, pattern="Lower", value = TRUE)]
#       ThetasLower <- c(-1e+21, ThetasLower[1] + c(0, cumsum(exp(ThetasLower[-1]))), 1e+21)
#       ThetasUpper <- beta[grep(par_thetas, pattern="Upper", value = TRUE)]
#       ThetasUpper <- c(-1e+21, ThetasUpper[1] + c(0, cumsum(exp(ThetasUpper[-1]))), 1e+21)
#       Thetas <- c(ThetasLower, ThetasUpper)
#   }
#   parammatrix[,c("th1", "th2")] <- cbind(Thetas[(DVs$response==1)*(nRatings+1) + DVs$rating],
#                                          Thetas[(DVs$response==1)*(nRatings+1) + DVs$rating+1])
#
#   # # Fill single missing muvis (index: 14) values with absolute value of decision drift v
#   # parammatrix[is.na(parammatrix[,14]), 14] <-
#   #   abs(parammatrix[is.na(parammatrix[,14]), 2])
#   #
#   # # Scale the parameters scaled by s by s
#   # # In the actual function:
#   # # cbind (a/s, v/s, t0, d, sz, sv/s, st0, z,
#   # #        tau, th1/s, th2/s, lambda, w, muvis/s, sigvis/s, svis/s, numeric_bounds)
#   # parammatrix[, c(1, 2, 7, 10, 11, 14, 15, 16)] <-
#   #   parammatrix[, c(1, 2,7, 10, 11, 14, 15, 16)]/parammatrix[,17]
#   #
#   # # bound between trial variability st0 to 2 (seconds)
#   # parammatrix[,8] <-  parammatrix[,8]*2
#   #
#   #
#   # # t0 <- t0+st0/2
#   # parammatrix[,3] <- parammatrix[,3] + parammatrix[,8]/2
#   #
#   # # probs <- d_WEVmu_fit(DVs$rt, DVs$response, parammatrix[,c(1:4, 6:8, 5, 9:16)], precision=2)
#   #
#   probs <- dWEV_parammatrix(DVs$rt, DVs$response, parammatrix, simult_conf,
#                             stop_on_error=TRUE, stop_on_zero=FALSE,
#                             precision = 3)
#   probs[probs==0] <- .Machine$double.xmin
#   logl <- sum(log(probs))
#   return(-logl)
#   # DVs <- cbind(DVs, parammatrix, probs)
# }
