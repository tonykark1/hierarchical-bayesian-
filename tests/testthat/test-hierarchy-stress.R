testthat::test_that("hierarchy scenarios span no signal to strong signal", {
  source("R/hierarchy_stress.R")

  s <- hierarchy_scenarios()
  testthat::expect_equal(s$scenario, c("none", "weak", "balanced", "strong"))
  testthat::expect_equal(s$sector_sd[1], 0)
  testthat::expect_gt(s$sector_sd[4], s$sector_sd[3])
  testthat::expect_lt(s$stock_sd[4], s$stock_sd[3])
})

testthat::test_that("sector corruption changes the requested stock-level mapping", {
  source("R/simulation.R")
  source("R/hierarchy_stress.R")

  sim <- simulate_factor_panel(n_stocks = 40, n_sectors = 5, n_obs = 26, seed = 31)
  corrupted <- corrupt_sector_labels(sim$panel, fraction = 0.25, seed = 32)
  acc <- sector_mapping_accuracy(corrupted)

  testthat::expect_equal(acc$n_stocks, 40)
  testthat::expect_equal(acc$incorrect, 10)
  testthat::expect_equal(acc$correct, 30)
  testthat::expect_equal(acc$accuracy, 0.75)

  mapping <- unique(corrupted[, c("stock", "sector")])
  testthat::expect_true(all(table(mapping$stock) == 1L))
})

testthat::test_that("zero corruption preserves all sector labels", {
  source("R/simulation.R")
  source("R/hierarchy_stress.R")

  sim <- simulate_factor_panel(n_stocks = 20, n_sectors = 4, n_obs = 26, seed = 33)
  corrupted <- corrupt_sector_labels(sim$panel, fraction = 0, seed = 34)
  acc <- sector_mapping_accuracy(corrupted)

  testthat::expect_equal(acc$accuracy, 1)
  testthat::expect_false(any(corrupted$sector_corrupted))
})

testthat::test_that("truth joins can ignore observed sector labels for corruption tests", {
  source("R/simulation.R")
  source("R/model_ols.R")
  source("R/hierarchy_stress.R")
  source("R/evaluation.R")

  sim <- simulate_factor_panel(n_stocks = 20, n_sectors = 4, n_obs = 52, seed = 35)
  corrupted <- corrupt_sector_labels(sim$panel, fraction = 0.5, seed = 36)
  est <- fit_ols_betas(corrupted, sim$metadata$factor_names)

  testthat::expect_error(
    join_estimates_to_truth(est, sim$truth$stock, match_sector = TRUE),
    "could not be matched"
  )

  joined <- join_estimates_to_truth(est, sim$truth$stock, match_sector = FALSE)
  testthat::expect_equal(nrow(joined), 20 * 6)
  testthat::expect_true(all(is.finite(joined$error)))
})

testthat::test_that("interval recovery reports coverage and width", {
  source("R/evaluation.R")

  d <- data.frame(
    factor = rep(c("MKT", "HML"), each = 2),
    beta = c(1.0, 1.2, 0.2, 0.4),
    lower = c(0.8, 1.0, 0.0, 0.5),
    upper = c(1.2, 1.4, 0.4, 0.7)
  )

  s <- summarise_interval_recovery(d)
  overall <- s[s$factor == "ALL", ]
  testthat::expect_equal(overall$coverage, 0.75)
  testthat::expect_true(overall$mean_width > 0)
})
