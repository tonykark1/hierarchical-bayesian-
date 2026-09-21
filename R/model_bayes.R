# Bayesian factor-model helpers.
#
# Phase 2 starts with the simplest defensible hierarchy:
#   market-level factor slope
#     + sector-level deviation
#     + stock-within-sector deviation
#
# Group-level correlations are disabled initially (|| in brms syntax) to avoid
# spending parameters on a large covariance structure before basic shrinkage is
# validated. Sector-specific variance parameters can be considered later if the
# simpler model is inadequate.

assert_brms_available <- function() {
  if (!requireNamespace("brms", quietly = TRUE)) {
    stop(
      "Package 'brms' is required for Bayesian model fitting. ",
      "Install brms and a Stan backend before calling this function."
    )
  }
  invisible(TRUE)
}

default_bayes_cores <- function(chains = 4L) {
  detected <- suppressWarnings(parallel::detectCores(logical = FALSE))
  if (length(detected) != 1L || is.na(detected) || detected < 1L) detected <- 1L
  max(1L, min(as.integer(chains), as.integer(detected)))
}

prepare_bayes_panel <- function(
    panel,
    factor_names = c("MKT", "SMB", "HML", "RMW", "CMA", "MOM")) {
  required <- c("stock", "sector", "excess_return", factor_names)
  missing_cols <- setdiff(required, names(panel))
  if (length(missing_cols) > 0L) {
    stop("panel is missing required columns: ", paste(missing_cols, collapse = ", "))
  }
  if (anyNA(panel[, required, drop = FALSE])) {
    stop("Bayesian model inputs must not contain missing values")
  }

  out <- panel
  if ("t" %in% names(out)) {
    out <- out[order(out$stock, out$t), , drop = FALSE]
  } else {
    out <- out[order(out$stock), , drop = FALSE]
  }

  out$stock <- factor(out$stock)
  out$sector <- factor(out$sector)
  out$sector_stock <- interaction(out$sector, out$stock, drop = TRUE, sep = "__")

  # A stock should map to one and only one sector in the v0.1 design.
  mapping <- unique(out[, c("stock", "sector"), drop = FALSE])
  counts <- table(mapping$stock)
  if (any(counts != 1L)) {
    stop("Each stock must belong to exactly one sector in the v0.1 hierarchy")
  }

  rownames(out) <- NULL
  out
}

make_hierarchical_factor_formula <- function(
    factor_names = c("MKT", "SMB", "HML", "RMW", "CMA", "MOM")) {
  if (length(factor_names) < 1L) stop("factor_names must contain at least one factor")
  if (anyDuplicated(factor_names)) stop("factor_names must be unique")

  fixed <- paste(factor_names, collapse = " + ")
  group_terms <- paste(c("1", factor_names), collapse = " + ")

  stats::as.formula(
    paste0(
      "excess_return ~ 1 + ", fixed,
      " + (", group_terms, " || sector)",
      " + (", group_terms, " || sector_stock)"
    )
  )
}

make_independent_factor_formula <- function(
    factor_names = c("MKT", "SMB", "HML", "RMW", "CMA", "MOM")) {
  if (length(factor_names) < 1L) stop("factor_names must contain at least one factor")
  stats::as.formula(
    paste("excess_return ~ 1 +", paste(factor_names, collapse = " + "))
  )
}

hierarchical_factor_priors <- function(
    factor_names = c("MKT", "SMB", "HML", "RMW", "CMA", "MOM"),
    beta_sd = 1.0,
    intercept_sd = 0.01,
    group_beta_sd = 0.35,
    group_intercept_sd = 0.01,
    residual_scale = 0.03) {
  assert_brms_available()

  stopifnot(
    beta_sd > 0,
    intercept_sd > 0,
    group_beta_sd > 0,
    group_intercept_sd > 0,
    residual_scale > 0
  )

  priors <- c(
    brms::set_prior(sprintf("normal(0, %s)", beta_sd), class = "b"),
    brms::set_prior(sprintf("normal(0, %s)", intercept_sd), class = "Intercept"),
    brms::set_prior(
      sprintf("student_t(3, 0, %s)", residual_scale),
      class = "sigma"
    ),
    brms::set_prior(
      sprintf("normal(0, %s)", group_intercept_sd),
      class = "sd",
      coef = "Intercept"
    )
  )

  for (f in factor_names) {
    priors <- c(
      priors,
      brms::set_prior(
        sprintf("normal(0, %s)", group_beta_sd),
        class = "sd",
        coef = f
      )
    )
  }

  priors
}

independent_factor_priors <- function(
    beta_sd = 1.0,
    intercept_sd = 0.01,
    residual_scale = 0.03) {
  assert_brms_available()
  stopifnot(beta_sd > 0, intercept_sd > 0, residual_scale > 0)

  c(
    brms::set_prior(sprintf("normal(0, %s)", beta_sd), class = "b"),
    brms::set_prior(sprintf("normal(0, %s)", intercept_sd), class = "Intercept"),
    brms::set_prior(
      sprintf("student_t(3, 0, %s)", residual_scale),
      class = "sigma"
    )
  )
}

fit_hierarchical_bayes <- function(
    panel,
    factor_names = c("MKT", "SMB", "HML", "RMW", "CMA", "MOM"),
    prior = NULL,
    chains = 4L,
    iter = 2000L,
    warmup = floor(iter / 2),
    cores = default_bayes_cores(chains),
    seed = 123L,
    backend = "cmdstanr",
    sample_prior = "yes",
    adapt_delta = 0.95,
    max_treedepth = 12L,
    refresh = 0L,
    ...) {
  assert_brms_available()
  d <- prepare_bayes_panel(panel, factor_names)
  formula <- make_hierarchical_factor_formula(factor_names)

  if (is.null(prior)) {
    prior <- hierarchical_factor_priors(factor_names)
  }

  brms::brm(
    formula = formula,
    data = d,
    family = brms::gaussian(),
    prior = prior,
    chains = chains,
    iter = iter,
    warmup = warmup,
    cores = cores,
    seed = seed,
    backend = backend,
    sample_prior = sample_prior,
    control = list(
      adapt_delta = adapt_delta,
      max_treedepth = max_treedepth
    ),
    refresh = refresh,
    ...
  )
}

fit_hierarchical_prior_predictive <- function(
    panel,
    factor_names = c("MKT", "SMB", "HML", "RMW", "CMA", "MOM"),
    prior = NULL,
    chains = 4L,
    iter = 1000L,
    warmup = floor(iter / 2),
    cores = default_bayes_cores(chains),
    seed = 123L,
    backend = "cmdstanr",
    adapt_delta = 0.95,
    max_treedepth = 12L,
    ...) {
  assert_brms_available()
  d <- prepare_bayes_panel(panel, factor_names)
  formula <- make_hierarchical_factor_formula(factor_names)

  if (is.null(prior)) {
    prior <- hierarchical_factor_priors(factor_names)
  }

  brms::brm(
    formula = formula,
    data = d,
    family = brms::gaussian(),
    prior = prior,
    chains = chains,
    iter = iter,
    warmup = warmup,
    cores = cores,
    seed = seed,
    backend = backend,
    sample_prior = "only",
    control = list(
      adapt_delta = adapt_delta,
      max_treedepth = max_treedepth
    ),
    refresh = 0L,
    ...
  )
}

# Independent Bayesian benchmark: one Bayesian regression per stock, with the
# same population-scale priors but no information sharing across stocks.
fit_independent_bayes <- function(
    panel,
    factor_names = c("MKT", "SMB", "HML", "RMW", "CMA", "MOM"),
    stocks = NULL,
    prior = NULL,
    chains = 2L,
    iter = 1500L,
    warmup = floor(iter / 2),
    cores = default_bayes_cores(chains),
    seed = 123L,
    backend = "cmdstanr",
    refresh = 0L,
    ...) {
  assert_brms_available()
  d <- prepare_bayes_panel(panel, factor_names)
  formula <- make_independent_factor_formula(factor_names)

  if (is.null(prior)) prior <- independent_factor_priors()
  if (is.null(stocks)) stocks <- levels(d$stock)
  stocks <- as.character(stocks)

  unknown <- setdiff(stocks, levels(d$stock))
  if (length(unknown) > 0L) {
    stop("Unknown stock(s): ", paste(unknown, collapse = ", "))
  }

  fits <- vector("list", length(stocks))
  names(fits) <- stocks

  for (i in seq_along(stocks)) {
    stock_i <- stocks[i]
    di <- d[as.character(d$stock) == stock_i, , drop = FALSE]

    fits[[i]] <- brms::brm(
      formula = formula,
      data = di,
      family = brms::gaussian(),
      prior = prior,
      chains = chains,
      iter = iter,
      warmup = warmup,
      cores = cores,
      seed = seed + i - 1L,
      backend = backend,
      sample_prior = "yes",
      refresh = refresh,
      ...
    )
  }

  structure(
    list(
      fits = fits,
      factor_names = factor_names,
      stocks = stocks
    ),
    class = "independent_bayes_factor_fits"
  )
}

summarise_draw_vector <- function(x, probs = c(0.025, 0.975)) {
  if (length(probs) != 2L || any(probs <= 0) || any(probs >= 1) || probs[1] >= probs[2]) {
    stop("probs must contain two ordered probabilities strictly between 0 and 1")
  }
  q <- stats::quantile(x, probs = probs, names = FALSE)
  c(
    estimate = mean(x),
    std_error = stats::sd(x),
    lower = q[1],
    upper = q[2]
  )
}

extract_independent_betas <- function(
    object,
    panel,
    factor_names = object$factor_names,
    probs = c(0.025, 0.975)) {
  assert_brms_available()
  if (!inherits(object, "independent_bayes_factor_fits")) {
    stop("object must come from fit_independent_bayes()")
  }

  d <- prepare_bayes_panel(panel, factor_names)
  rows <- vector("list", length(object$stocks))

  for (i in seq_along(object$stocks)) {
    stock_i <- object$stocks[i]
    fit <- object$fits[[stock_i]]
    draws <- brms::fixef(fit, summary = FALSE)
    sector_i <- as.character(d$sector[as.character(d$stock) == stock_i][1])
    n_i <- sum(as.character(d$stock) == stock_i)

    stock_rows <- lapply(factor_names, function(f) {
      if (!f %in% colnames(draws)) stop("Factor not found in posterior draws: ", f)
      s <- summarise_draw_vector(draws[, f], probs)
      data.frame(
        stock = stock_i,
        sector = sector_i,
        factor = f,
        estimate = unname(s["estimate"]),
        std_error = unname(s["std_error"]),
        lower = unname(s["lower"]),
        upper = unname(s["upper"]),
        n_obs = n_i,
        estimator = "Independent Bayes",
        stringsAsFactors = FALSE
      )
    })
    rows[[i]] <- do.call(rbind, stock_rows)
  }

  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

# Reconstruct stock-level slope draws as:
# fixed market slope + sector deviation + stock-within-sector deviation.
extract_hierarchical_betas <- function(
    fit,
    panel,
    factor_names = c("MKT", "SMB", "HML", "RMW", "CMA", "MOM"),
    probs = c(0.025, 0.975)) {
  assert_brms_available()
  d <- prepare_bayes_panel(panel, factor_names)
  mapping <- unique(d[, c("stock", "sector", "sector_stock"), drop = FALSE])
  mapping$stock <- as.character(mapping$stock)
  mapping$sector <- as.character(mapping$sector)
  mapping$sector_stock <- as.character(mapping$sector_stock)

  fixed <- brms::fixef(fit, summary = FALSE)
  random <- brms::ranef(fit, summary = FALSE)

  if (!all(c("sector", "sector_stock") %in% names(random))) {
    stop("Expected sector and sector_stock random effects were not found")
  }

  sector_re <- random[["sector"]]
  stock_re <- random[["sector_stock"]]

  sector_levels <- dimnames(sector_re)[[2]]
  sector_coefs <- dimnames(sector_re)[[3]]
  stock_levels <- dimnames(stock_re)[[2]]
  stock_coefs <- dimnames(stock_re)[[3]]

  rows <- vector("list", nrow(mapping))

  for (i in seq_len(nrow(mapping))) {
    stock_i <- mapping$stock[i]
    sector_i <- mapping$sector[i]
    nested_i <- mapping$sector_stock[i]
    n_i <- sum(as.character(d$stock) == stock_i)

    if (!sector_i %in% sector_levels) stop("Sector level missing from ranef draws: ", sector_i)
    if (!nested_i %in% stock_levels) stop("Stock level missing from ranef draws: ", nested_i)

    stock_rows <- lapply(factor_names, function(f) {
      if (!f %in% colnames(fixed)) stop("Factor missing from fixed-effect draws: ", f)
      if (!f %in% sector_coefs) stop("Factor missing from sector random effects: ", f)
      if (!f %in% stock_coefs) stop("Factor missing from stock random effects: ", f)

      beta_draws <- fixed[, f] +
        sector_re[, sector_i, f] +
        stock_re[, nested_i, f]
      s <- summarise_draw_vector(beta_draws, probs)

      data.frame(
        stock = stock_i,
        sector = sector_i,
        factor = f,
        estimate = unname(s["estimate"]),
        std_error = unname(s["std_error"]),
        lower = unname(s["lower"]),
        upper = unname(s["upper"]),
        n_obs = n_i,
        estimator = "Hierarchical Bayes",
        stringsAsFactors = FALSE
      )
    })

    rows[[i]] <- do.call(rbind, stock_rows)
  }

  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

summarise_brms_diagnostics <- function(fit, max_treedepth = 12L) {
  assert_brms_available()
  if (!requireNamespace("posterior", quietly = TRUE)) {
    stop("Package 'posterior' is required for diagnostic summaries")
  }

  draws <- posterior::as_draws_array(fit)
  diag <- posterior::summarise_draws(draws)
  nuts <- brms::nuts_params(fit)

  rhat_values <- diag$rhat[is.finite(diag$rhat)]
  bulk_values <- diag$ess_bulk[is.finite(diag$ess_bulk)]
  tail_values <- diag$ess_tail[is.finite(diag$ess_tail)]

  divergences <- sum(
    nuts$Parameter == "divergent__" & nuts$Value > 0,
    na.rm = TRUE
  )
  treedepth_hits <- sum(
    nuts$Parameter == "treedepth__" & nuts$Value >= max_treedepth,
    na.rm = TRUE
  )

  data.frame(
    max_rhat = if (length(rhat_values)) max(rhat_values) else NA_real_,
    min_bulk_ess = if (length(bulk_values)) min(bulk_values) else NA_real_,
    min_tail_ess = if (length(tail_values)) min(tail_values) else NA_real_,
    divergences = divergences,
    max_treedepth_hits = treedepth_hits,
    stringsAsFactors = FALSE
  )
}

# Lightweight structural checks that do not require fitting a Stan model.
validate_hierarchical_design <- function(
    panel,
    factor_names = c("MKT", "SMB", "HML", "RMW", "CMA", "MOM")) {
  d <- prepare_bayes_panel(panel, factor_names)

  n_stocks <- nlevels(d$stock)
  n_sectors <- nlevels(d$sector)
  mapping <- unique(d[, c("stock", "sector"), drop = FALSE])
  stocks_per_sector <- table(mapping$sector)

  data.frame(
    n_rows = nrow(d),
    n_stocks = n_stocks,
    n_sectors = n_sectors,
    min_stocks_per_sector = min(stocks_per_sector),
    median_stocks_per_sector = stats::median(as.numeric(stocks_per_sector)),
    max_stocks_per_sector = max(stocks_per_sector),
    n_factors = length(factor_names),
    stringsAsFactors = FALSE
  )
}
