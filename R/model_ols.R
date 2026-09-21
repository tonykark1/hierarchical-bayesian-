# OLS baseline for stock-level factor exposures.

fit_ols_betas <- function(panel,
                          factor_names = c("MKT", "SMB", "HML", "RMW", "CMA", "MOM")) {
  required <- c("stock", "sector", "excess_return", factor_names)
  missing_cols <- setdiff(required, names(panel))
  if (length(missing_cols) > 0L) {
    stop("panel is missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  stocks <- unique(panel$stock)
  out <- vector("list", length(stocks))

  rhs <- paste(factor_names, collapse = " + ")
  fml <- stats::as.formula(paste("excess_return ~", rhs))

  for (i in seq_along(stocks)) {
    d <- panel[panel$stock == stocks[i], , drop = FALSE]
    fit <- stats::lm(fml, data = d)
    sm <- summary(fit)$coefficients

    beta_rows <- sm[factor_names, , drop = FALSE]
    out[[i]] <- data.frame(
      stock = stocks[i],
      sector = d$sector[1],
      factor = factor_names,
      estimate = beta_rows[, "Estimate"],
      std_error = beta_rows[, "Std. Error"],
      t_value = beta_rows[, "t value"],
      n_obs = stats::nobs(fit),
      stringsAsFactors = FALSE
    )
  }

  ans <- do.call(rbind, out)
  rownames(ans) <- NULL
  ans
}
