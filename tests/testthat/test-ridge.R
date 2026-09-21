testthat::test_that("expanding splits never use future observations in training", {
  source("R/model_ridge.R")

  splits <- make_expanding_splits(n = 52, initial = 26, assess = 5, max_splits = 4)
  testthat::expect_true(length(splits) >= 1)

  for (s in splits) {
    testthat::expect_lt(max(s$analysis), min(s$assessment))
    testthat::expect_equal(s$analysis, seq_len(max(s$analysis)))
    testthat::expect_false(any(s$analysis %in% s$assessment))
  }
})

testthat::test_that("ridge with near-zero penalty approximates OLS", {
  source("R/simulation.R")
  source("R/model_ols.R")
  source("R/model_ridge.R")

  sim <- simulate_factor_panel(n_stocks = 8, n_sectors = 2, n_obs = 156, seed = 11)
  stock <- sim$panel$stock[1]
  d <- sim$panel[sim$panel$stock == stock, , drop = FALSE]
  factors <- sim$metadata$factor_names

  ols <- fit_ols_betas(d, factors)
  ridge <- fit_ridge_core(
    as.matrix(d[, factors, drop = FALSE]),
    d$excess_return,
    lambda = 1e-10
  )

  testthat::expect_equal(
    unname(ridge$beta[factors]),
    ols$estimate,
    tolerance = 1e-6
  )
})

testthat::test_that("ridge tuning returns a lambda from the candidate grid", {
  source("R/simulation.R")
  source("R/model_ridge.R")

  sim <- simulate_factor_panel(n_stocks = 8, n_sectors = 2, n_obs = 52, seed = 12)
  d <- sim$panel[sim$panel$stock == sim$panel$stock[1], , drop = FALSE]
  grid <- c(1e-3, 1e-2, 1e-1, 1, 10)

  tuned <- select_ridge_lambda(
    d,
    factor_names = sim$metadata$factor_names,
    lambda_grid = grid,
    initial = 26,
    assess = 5,
    max_splits = 3
  )

  testthat::expect_true(tuned$lambda %in% grid)
  testthat::expect_true(is.finite(tuned$cv_mse))
  testthat::expect_equal(nrow(tuned$grid), length(grid))
})

testthat::test_that("ridge estimator has OLS-compatible beta schema and truth join", {
  source("R/simulation.R")
  source("R/model_ridge.R")
  source("R/evaluation.R")

  sim <- simulate_factor_panel(n_stocks = 12, n_sectors = 3, n_obs = 52, seed = 13)
  est <- fit_ridge_betas(
    sim$panel,
    sim$metadata$factor_names,
    lambda_grid = c(1e-3, 1e-2, 1e-1, 1, 10),
    max_splits = 3
  )

  expected_core <- c("stock", "sector", "factor", "estimate", "std_error", "t_value", "n_obs")
  testthat::expect_true(all(expected_core %in% names(est)))
  testthat::expect_equal(nrow(est), 12 * 6)
  testthat::expect_true(all(est$lambda > 0))

  joined <- join_estimates_to_truth(est, sim$truth$stock)
  testthat::expect_equal(nrow(joined), 12 * 6)
  testthat::expect_false(anyNA(joined$beta))
  testthat::expect_true(all(is.finite(joined$error)))
})

testthat::test_that("ridge rejects invalid predictor scales and lambdas", {
  source("R/model_ridge.R")

  x <- cbind(a = rep(1, 10), b = seq_len(10))
  y <- rnorm(10)
  testthat::expect_error(fit_ridge_core(x, y, 1), "positive finite")
  testthat::expect_error(fit_ridge_core(cbind(a = 1:10, b = 11:20), y, -1), "lambda")
})
