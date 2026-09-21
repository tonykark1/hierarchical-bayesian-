# Evaluation helpers shared by simulation and later out-of-sample experiments.

join_estimates_to_truth <- function(estimates, truth) {
  req_est <- c("stock", "sector", "factor", "estimate")
  req_truth <- c("stock", "sector", "factor", "beta")

  miss_est <- setdiff(req_est, names(estimates))
  miss_truth <- setdiff(req_truth, names(truth))
  if (length(miss_est) > 0L) {
    stop("estimates missing: ", paste(miss_est, collapse = ", "))
  }
  if (length(miss_truth) > 0L) {
    stop("truth missing: ", paste(miss_truth, collapse = ", "))
  }

  merged <- merge(
    estimates,
    truth,
    by = c("stock", "sector", "factor"),
    all.x = TRUE,
    sort = FALSE
  )

  if (anyNA(merged$beta)) {
    stop("Some estimates could not be matched to true betas")
  }

  merged$error <- merged$estimate - merged$beta
  merged$abs_error <- abs(merged$error)
  merged$sq_error <- merged$error^2
  merged
}

summarise_beta_recovery <- function(joined) {
  required <- c("factor", "error", "abs_error", "sq_error")
  missing_cols <- setdiff(required, names(joined))
  if (length(missing_cols) > 0L) {
    stop("joined data missing: ", paste(missing_cols, collapse = ", "))
  }

  factors <- unique(joined$factor)
  rows <- lapply(factors, function(f) {
    d <- joined[joined$factor == f, , drop = FALSE]
    data.frame(
      factor = f,
      n = nrow(d),
      bias = mean(d$error),
      mae = mean(d$abs_error),
      rmse = sqrt(mean(d$sq_error)),
      stringsAsFactors = FALSE
    )
  })

  overall <- data.frame(
    factor = "ALL",
    n = nrow(joined),
    bias = mean(joined$error),
    mae = mean(joined$abs_error),
    rmse = sqrt(mean(joined$sq_error)),
    stringsAsFactors = FALSE
  )

  out <- rbind(do.call(rbind, rows), overall)
  rownames(out) <- NULL
  out
}

run_ols_recovery_experiment <- function(
    sample_sizes = c(26L, 52L, 104L, 156L, 260L, 520L),
    n_stocks = 100L,
    n_sectors = 10L,
    seed = 123L) {

  results <- vector("list", length(sample_sizes))

  for (j in seq_along(sample_sizes)) {
    sim <- simulate_factor_panel(
      n_stocks = n_stocks,
      n_sectors = n_sectors,
      n_obs = sample_sizes[j],
      seed = seed + j - 1L
    )

    est <- fit_ols_betas(sim$panel, sim$metadata$factor_names)
    joined <- join_estimates_to_truth(est, sim$truth$stock)
    summary <- summarise_beta_recovery(joined)
    summary$n_obs <- sample_sizes[j]
    summary$estimator <- "OLS"
    results[[j]] <- summary
  }

  out <- do.call(rbind, results)
  rownames(out) <- NULL
  out
}

run_ridge_recovery_experiment <- function(
    sample_sizes = c(26L, 52L, 104L, 156L, 260L, 520L),
    n_stocks = 100L,
    n_sectors = 10L,
    seed = 123L,
    lambda_grid = ridge_lambda_grid()) {

  results <- vector("list", length(sample_sizes))

  for (j in seq_along(sample_sizes)) {
    sim <- simulate_factor_panel(
      n_stocks = n_stocks,
      n_sectors = n_sectors,
      n_obs = sample_sizes[j],
      seed = seed + j - 1L
    )

    est <- fit_ridge_betas(
      sim$panel,
      sim$metadata$factor_names,
      lambda_grid = lambda_grid
    )
    joined <- join_estimates_to_truth(est, sim$truth$stock)
    summary <- summarise_beta_recovery(joined)
    summary$n_obs <- sample_sizes[j]
    summary$estimator <- "Ridge"
    results[[j]] <- summary
  }

  out <- do.call(rbind, results)
  rownames(out) <- NULL
  out
}

run_baseline_recovery_experiment <- function(
    sample_sizes = c(26L, 52L, 104L, 156L, 260L, 520L),
    n_stocks = 100L,
    n_sectors = 10L,
    seed = 123L,
    lambda_grid = ridge_lambda_grid()) {

  ols <- run_ols_recovery_experiment(
    sample_sizes = sample_sizes,
    n_stocks = n_stocks,
    n_sectors = n_sectors,
    seed = seed
  )

  ridge <- run_ridge_recovery_experiment(
    sample_sizes = sample_sizes,
    n_stocks = n_stocks,
    n_sectors = n_sectors,
    seed = seed,
    lambda_grid = lambda_grid
  )

  out <- rbind(ols, ridge)
  out <- out[, c("estimator", "n_obs", "factor", "n", "bias", "mae", "rmse")]
  rownames(out) <- NULL
  out
}
