testthat::test_that("simulation is reproducible under a fixed seed", {
  source("R/simulation.R")

  a <- simulate_factor_panel(n_stocks = 20, n_sectors = 4, n_obs = 52, seed = 42)
  b <- simulate_factor_panel(n_stocks = 20, n_sectors = 4, n_obs = 52, seed = 42)

  testthat::expect_identical(a$panel, b$panel)
  testthat::expect_identical(a$truth$stock, b$truth$stock)
  testthat::expect_identical(a$truth$sector, b$truth$sector)
})

testthat::test_that("simulation dimensions and hierarchy are correct", {
  source("R/simulation.R")

  n_stocks <- 24
  n_sectors <- 6
  n_obs <- 104
  k <- 6

  sim <- simulate_factor_panel(
    n_stocks = n_stocks,
    n_sectors = n_sectors,
    n_obs = n_obs,
    seed = 7
  )

  testthat::expect_equal(nrow(sim$panel), n_stocks * n_obs)
  testthat::expect_equal(length(unique(sim$panel$stock)), n_stocks)
  testthat::expect_equal(length(unique(sim$panel$sector)), n_sectors)
  testthat::expect_equal(nrow(sim$truth$stock), n_stocks * k)
  testthat::expect_equal(nrow(sim$truth$sector), n_sectors * k)
})

testthat::test_that("OLS estimates join exactly to true betas", {
  source("R/simulation.R")
  source("R/model_ols.R")
  source("R/evaluation.R")

  sim <- simulate_factor_panel(n_stocks = 20, n_sectors = 4, n_obs = 104, seed = 9)
  est <- fit_ols_betas(sim$panel, sim$metadata$factor_names)
  joined <- join_estimates_to_truth(est, sim$truth$stock)

  testthat::expect_equal(nrow(joined), 20 * 6)
  testthat::expect_false(anyNA(joined$beta))
  testthat::expect_true(all(is.finite(joined$error)))
})
