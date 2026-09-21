testthat::test_that("Bayesian panel preparation preserves one sector per stock", {
  source("R/simulation.R")
  source("R/model_bayes.R")

  sim <- simulate_factor_panel(n_stocks = 24, n_sectors = 6, n_obs = 52, seed = 21)
  d <- prepare_bayes_panel(sim$panel, sim$metadata$factor_names)

  testthat::expect_true(is.factor(d$stock))
  testthat::expect_true(is.factor(d$sector))
  testthat::expect_true(is.factor(d$sector_stock))
  testthat::expect_equal(nlevels(d$stock), 24)
  testthat::expect_equal(nlevels(d$sector), 6)
  testthat::expect_equal(nlevels(d$sector_stock), 24)

  mapping <- unique(d[, c("stock", "sector")])
  testthat::expect_true(all(table(mapping$stock) == 1L))
})

testthat::test_that("hierarchical formula contains market, sector and stock levels", {
  source("R/model_bayes.R")

  factors <- c("MKT", "HML", "MOM")
  txt <- paste(deparse(make_hierarchical_factor_formula(factors)), collapse = " ")

  testthat::expect_match(txt, "excess_return")
  testthat::expect_match(txt, "MKT")
  testthat::expect_match(txt, "HML")
  testthat::expect_match(txt, "MOM")
  testthat::expect_match(txt, "sector")
  testthat::expect_match(txt, "sector_stock")
  testthat::expect_match(txt, "\\|\\|")
})

testthat::test_that("independent formula has no pooling terms", {
  source("R/model_bayes.R")

  factors <- c("MKT", "HML", "MOM")
  txt <- paste(deparse(make_independent_factor_formula(factors)), collapse = " ")

  testthat::expect_match(txt, "MKT")
  testthat::expect_match(txt, "HML")
  testthat::expect_match(txt, "MOM")
  testthat::expect_false(grepl("sector", txt, fixed = TRUE))
  testthat::expect_false(grepl("|", txt, fixed = TRUE))
})

testthat::test_that("Bayesian panel preparation rejects changing sector membership", {
  source("R/simulation.R")
  source("R/model_bayes.R")

  sim <- simulate_factor_panel(n_stocks = 12, n_sectors = 3, n_obs = 26, seed = 22)
  bad <- sim$panel
  first_stock <- bad$stock[1]
  idx <- which(bad$stock == first_stock)[1]
  bad$sector[idx] <- "OTHER"

  testthat::expect_error(
    prepare_bayes_panel(bad, sim$metadata$factor_names),
    "exactly one sector"
  )
})

testthat::test_that("hierarchical design summary reports expected dimensions", {
  source("R/simulation.R")
  source("R/model_bayes.R")

  sim <- simulate_factor_panel(n_stocks = 30, n_sectors = 5, n_obs = 26, seed = 23)
  summary <- validate_hierarchical_design(sim$panel, sim$metadata$factor_names)

  testthat::expect_equal(summary$n_stocks, 30)
  testthat::expect_equal(summary$n_sectors, 5)
  testthat::expect_equal(summary$n_factors, 6)
  testthat::expect_equal(summary$min_stocks_per_sector, 6)
  testthat::expect_equal(summary$max_stocks_per_sector, 6)
})

testthat::test_that("draw-vector summary returns posterior-compatible schema", {
  source("R/model_bayes.R")

  x <- seq(-1, 1, length.out = 101)
  s <- summarise_draw_vector(x)

  testthat::expect_equal(names(s), c("estimate", "std_error", "lower", "upper"))
  testthat::expect_equal(unname(s["estimate"]), 0, tolerance = 1e-12)
  testthat::expect_lt(s["lower"], s["estimate"])
  testthat::expect_gt(s["upper"], s["estimate"])
})

testthat::test_that("default core selection is always a positive integer", {
  source("R/model_bayes.R")

  cores <- default_bayes_cores(4L)
  testthat::expect_true(is.numeric(cores))
  testthat::expect_gte(cores, 1L)
  testthat::expect_lte(cores, 4L)
})
