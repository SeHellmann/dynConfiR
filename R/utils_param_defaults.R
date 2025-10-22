#' Fill optional parameter columns with default values.
#'
#' Ensures that optional scalar parameters exist and replaces `NA` entries with
#' the provided defaults. This keeps the original column order intact (new
#' columns are appended) and works for one-row data frames and list inputs that
#' can be coerced to a data frame.
#'
#' @param paramDf Data frame or list specifying model parameters.
#' @param defaults Named list or vector of default values.
#'
#' @return A data frame with missing/`NA` optional columns replaced.
#' @keywords internal
fill_optional_params <- function(paramDf, defaults) {
  paramDf <- as.data.frame(paramDf, stringsAsFactors = FALSE)
  if (length(defaults) == 0L || nrow(paramDf) == 0L) {
    return(paramDf)
  }

  for (nm in names(defaults)) {
    default_val <- defaults[[nm]]
    if (!nm %in% names(paramDf)) {
      paramDf[[nm]] <- rep_len(default_val, nrow(paramDf))
      next
    }

    missing_idx <- is.na(paramDf[[nm]])
    if (any(missing_idx)) {
      paramDf[[nm]][missing_idx] <- rep_len(default_val, sum(missing_idx))
    }
  }

  paramDf
}
