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

fit_hierarchical_bayes <- function(
    panel,
    factor_names = c("MKT", "SMB", "HML", "RMW", "CMA", "MOM"),
    prior = NULL,
    chains = 4L,
    iter = 2000L,
    warmup = floor(iter / 2),
    cores = min(chains, parallel::detectCores(logical = FALSE)),
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
    seed = 123L,
    backend = "cmdstanr",
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
    warmup = 0L,
    seed = seed,
    backend = backend,
    sample_prior = "only",
    refresh = 0L,
    ...
  )
}

# Lightweight structural checks that do not require fitting a Stan model.
validate_hierarchical_design <- function(
    panel,
    factor_names = c("MKT", "SMB", "HML", "RMW", "CMA", "MOM")) {
  d <- prepare_bayes_panel(panel, factor_names)

  n_stocks <- nlevels(d$stock)
  n_sectors <- nlevels(d$sector)
  stocks_per_sector <- table(unique(d[, c("stock", "sector")])$sector)

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
