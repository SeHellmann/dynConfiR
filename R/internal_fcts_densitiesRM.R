### ToDo: Maybe put all parameters in arguments (as in dWEV etc.)
### Alternatively use params as argument for prepare_RaceModel_parameter
prepare_RaceModel_parameter <- function(response,mu1, mu2,
                                        a, b, s, th1, th2,t0, st0, wx, wrt, wint, nn) {

  if ( (length(mu1) == 1) &
       (length(mu2) == 1) &
       (length(a) == 1) &
       (length(b) == 1) &
       (length(s) == 1) &
       (length(th1) == 1) &
       (length(th2) == 1) &
       (length(t0) == 1) &
       (length(st0) == 1) &
       (length(wx) == 1) &
       (length(wrt) == 1) &
       (length(wint) == 1)) {
    skip_checks <- TRUE
  } else {
    skip_checks <- FALSE
  }

  response <- as.numeric(response)
  if (any(!(response %in% 1:2)))
    stop("response needs to be  %in% 1:2!")
  numeric_bounds <- as.integer(response)

  numeric_bounds <- rep(numeric_bounds, length.out = nn)
  if (!skip_checks) {
    # all parameters brought to length of rt
    mu1 <- rep(mu1, length.out = nn)
    mu2 <- rep(mu2, length.out = nn)
    b <- rep(b, length.out = nn)
    s <- rep(s, length.out = nn)
    a <- rep(a, length.out = nn)
    th1 <- rep(th1, length.out = nn)
    th2 <- rep(th2, length.out = nn)
    t0 <- rep(t0, length.out = nn)
    st0 <- rep(st0, length.out = nn)
    wx <- rep(wx, length.out = nn)
    wrt <- rep(wrt, length.out = nn)
    wint <- rep(wint, length.out = nn)
  }
  th1[th1==-Inf] <- 0
  th2[th2==Inf] <- .Machine$double.xmax
  # Build parameter matrix (and divide a, v, and sv, by s)
  params <- cbind (mu1, mu2, -a, -b, s, th1, th2, st0, wx, wrt, wint, t0, numeric_bounds)

  # Check for illegal parameter values
  if(ncol(params)<13) stop("Not enough parameters supplied: probable attempt to pass NULL values?")
  if(!is.numeric(params)) stop("Parameters need to be numeric.")
  if (any(is.na(params)) || !all(is.finite(params))) {
    stop("Parameters need to be numeric and finite.")
  }


  if (!skip_checks) {
    parameter_char <- apply(params, 1, paste0, collapse = "\t")
    parameter_factor <- factor(parameter_char, levels = unique(parameter_char))
    parameter_indices <- split(seq_len(nn), f = parameter_factor)
  } else {
    if (all(numeric_bounds == 2L) | all(numeric_bounds == 1L)) {
      parameter_indices <- list(
        seq_len(nn)
      )
    } else {
      parameter_indices <- list(
        seq_len(nn)[numeric_bounds == 2L],
        seq_len(nn)[numeric_bounds == 1L]
      )
    }
  }
  list(
    params = params
    , parameter_indices = parameter_indices
  )
}
