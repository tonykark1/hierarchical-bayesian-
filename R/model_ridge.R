# Leakage-safe ridge benchmark for stock-level factor exposures.
#
# This implementation intentionally uses base R only. Predictors are standardised
# using each training fold only, the intercept is never penalised, and lambda is
# chosen with expanding-window validation so future observations cannot leak into
# the tuning sample.

ridge_lambda_grid <- function(n = 31L, min_lambda = 1e-4, max_lambda = 1e3) {
  stopifnot(n >= 2L, min_lambda > 0, max_lambda > min_lambda)
  exp(seq(log(min_lambda), log(max_lambda), length.out = n))
}

make_expanding_splits <- function(n,
                                  initial = max(12L, floor(0.5 * n)),
                                  assess = max(1L, floor(0.1 * n)),
                                  max_splits = 5L) {
  n <- as.integer(n)
  initial <- as.integer(initial)
  assess <- as.integer(assess)
  max_splits <- as.integer(max_splits)

  if (n < 4L) stop("n must be at least 4")
  if (initial < 2L || initial >= n) stop("initial must be >= 2 and < n")
  if (assess < 1L) stop("assess must be >= 1")
  if (initial + assess > n) stop("initial + assess must not exceed n")
  if (max_splits < 1L) stop("max_splits must be >= 1")

  train_ends <- seq.int(initial, n - assess, by = assess)
  if (length(train_ends) > max_splits) {
    train_ends <- tail(train_ends, max_splits)
  }

  lapply(train_ends, function(train_end) {
    list(
      analysis = seq_len(train_end),
      assessment = seq.int(train_end + 1L, train_end + assess)
    )
  })
}

fit_ridge_core <- function(x, y, lambda) {
  x <- as.matrix(x)
  y <- as.numeric(y)
  lambda <- as.numeric(lambda)

  if (nrow(x) != length(y)) stop("x and y have incompatible dimensions")
  if (nrow(x) < 2L) stop("Need at least two observations")
  if (ncol(x) < 1L) stop("Need at least one predictor")
  if (!is.finite(lambda) || lambda < 0) stop("lambda must be finite and non-negative")
  if (anyNA(x) || anyNA(y)) stop("x and y must not contain missing values")

  x_mean <- colMeans(x)
  x_sd <- apply(x, 2L, stats::sd)
  if (any(!is.finite(x_sd)) || any(x_sd <= 0)) {
    stop("Every predictor must have positive finite training-sample standard deviation")
  }

  y_mean <- mean(y)
  x_std <- sweep(sweep(x, 2L, x_mean, FUN = "-"), 2L, x_sd, FUN = "/")
  y_ctr <- y - y_mean

  penalty <- diag(lambda, nrow = ncol(x_std), ncol = ncol(x_std))
  beta_std <- drop(solve(crossprod(x_std) + penalty, crossprod(x_std, y_ctr)))
  beta <- beta_std / x_sd
  intercept <- y_mean - sum(beta * x_mean)

  list(
    intercept = unname(intercept),
    beta = stats::setNames(unname(beta), colnames(x)),
    x_mean = x_mean,
    x_sd = x_sd,
    lambda = lambda
  )
}

predict_ridge_core <- function(fit, x) {
  x <- as.matrix(x)
  if (!all(names(fit$beta) %in% colnames(x))) {
    stop("x is missing predictors used by the fitted ridge model")
  }
  x <- x[, names(fit$beta), drop = FALSE]
  drop(fit$intercept + x %*% fit$beta)
}

select_ridge_lambda <- function(d,
                                factor_names,
                                lambda_grid = ridge_lambda_grid(),
                                initial = max(12L, floor(0.5 * nrow(d))),
                                assess = max(1L, floor(0.1 * nrow(d))),
                                max_splits = 5L) {
  required <- c("excess_return", factor_names)
  missing_cols <- setdiff(required, names(d))
  if (length(missing_cols) > 0L) {
    stop("d is missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  if ("t" %in% names(d)) d <- d[order(d$t), , drop = FALSE]
  splits <- make_expanding_splits(
    n = nrow(d), initial = initial, assess = assess, max_splits = max_splits
  )

  losses <- matrix(NA_real_, nrow = length(lambda_grid), ncol = length(splits))

  for (s in seq_along(splits)) {
    train <- d[splits[[s]]$analysis, , drop = FALSE]
    valid <- d[splits[[s]]$assessment, , drop = FALSE]

    x_train <- as.matrix(train[, factor_names, drop = FALSE])
    y_train <- train$excess_return
    x_valid <- as.matrix(valid[, factor_names, drop = FALSE])
    y_valid <- valid$excess_return

    for (j in seq_along(lambda_grid)) {
      fit <- fit_ridge_core(x_train, y_train, lambda_grid[j])
      pred <- predict_ridge_core(fit, x_valid)
      losses[j, s] <- mean((y_valid - pred)^2)
    }
  }

  cv_mse <- rowMeans(losses)
  best <- which.min(cv_mse)

  list(
    lambda = lambda_grid[best],
    cv_mse = cv_mse[best],
    grid = data.frame(lambda = lambda_grid, cv_mse = cv_mse),
    splits = splits
  )
}

fit_ridge_betas <- function(panel,
                            factor_names = c("MKT", "SMB", "HML", "RMW", "CMA", "MOM"),
                            lambda_grid = ridge_lambda_grid(),
                            initial_fraction = 0.5,
                            assess_fraction = 0.1,
                            max_splits = 5L) {
  required <- c("stock", "sector", "excess_return", factor_names)
  missing_cols <- setdiff(required, names(panel))
  if (length(missing_cols) > 0L) {
    stop("panel is missing required columns: ", paste(missing_cols, collapse = ", "))
  }
  if (initial_fraction <= 0 || initial_fraction >= 1) {
    stop("initial_fraction must be in (0, 1)")
  }
  if (assess_fraction <= 0 || assess_fraction >= 1) {
    stop("assess_fraction must be in (0, 1)")
  }

  stocks <- unique(panel$stock)
  out <- vector("list", length(stocks))

  for (i in seq_along(stocks)) {
    d <- panel[panel$stock == stocks[i], , drop = FALSE]
    if ("t" %in% names(d)) d <- d[order(d$t), , drop = FALSE]

    n <- nrow(d)
    initial <- max(length(factor_names) + 3L, floor(initial_fraction * n))
    assess <- max(1L, floor(assess_fraction * n))
    if (initial + assess > n) assess <- n - initial
    if (assess < 1L) stop("Not enough observations for ridge validation")

    tuning <- select_ridge_lambda(
      d = d,
      factor_names = factor_names,
      lambda_grid = lambda_grid,
      initial = initial,
      assess = assess,
      max_splits = max_splits
    )

    fit <- fit_ridge_core(
      x = as.matrix(d[, factor_names, drop = FALSE]),
      y = d$excess_return,
      lambda = tuning$lambda
    )

    out[[i]] <- data.frame(
      stock = stocks[i],
      sector = d$sector[1],
      factor = factor_names,
      estimate = unname(fit$beta[factor_names]),
      std_error = NA_real_,
      t_value = NA_real_,
      n_obs = n,
      lambda = tuning$lambda,
      cv_mse = tuning$cv_mse,
      stringsAsFactors = FALSE
    )
  }

  ans <- do.call(rbind, out)
  rownames(ans) <- NULL
  ans
}
