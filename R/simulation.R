# Synthetic data generator for the hierarchical factor lab.
#
# The first phase intentionally uses base R only. This keeps the DGP transparent
# and makes parameter-recovery failures easier to diagnose before adding brms or
# cmdstanr.

simulate_factor_panel <- function(
    n_stocks = 100L,
    n_sectors = 10L,
    n_obs = 156L,
    factor_names = c("MKT", "SMB", "HML", "RMW", "CMA", "MOM"),
    market_beta_mean = c(MKT = 1.00, SMB = 0.20, HML = 0.20,
                         RMW = 0.10, CMA = 0.10, MOM = 0.10),
    sector_sd = 0.20,
    stock_sd = 0.15,
    residual_sd_range = c(0.015, 0.035),
    factor_sd = c(MKT = 0.025, SMB = 0.015, HML = 0.015,
                  RMW = 0.012, CMA = 0.012, MOM = 0.018),
    factor_correlation = 0.10,
    seed = 123L) {

  stopifnot(
    n_stocks >= 2L,
    n_sectors >= 1L,
    n_sectors <= n_stocks,
    n_obs >= 10L,
    length(residual_sd_range) == 2L,
    residual_sd_range[1] > 0,
    residual_sd_range[2] >= residual_sd_range[1],
    factor_correlation > -1,
    factor_correlation < 1
  )

  if (!all(factor_names %in% names(market_beta_mean))) {
    stop("market_beta_mean must contain every factor in factor_names")
  }
  if (!all(factor_names %in% names(factor_sd))) {
    stop("factor_sd must contain every factor in factor_names")
  }

  set.seed(seed)

  k <- length(factor_names)
  stock_id <- sprintf("STK%03d", seq_len(n_stocks))
  sector_id <- sprintf("SEC%02d", seq_len(n_sectors))
  stock_sector <- sector_id[rep(seq_len(n_sectors), length.out = n_stocks)]

  market_mean <- market_beta_mean[factor_names]

  # Sector betas are draws around the market-level factor exposure.
  sector_beta_matrix <- matrix(
    rep(market_mean, each = n_sectors),
    nrow = n_sectors,
    ncol = k,
    dimnames = list(sector_id, factor_names)
  ) + matrix(
    rnorm(n_sectors * k, mean = 0, sd = sector_sd),
    nrow = n_sectors,
    ncol = k
  )

  # Stock betas are draws around their sector-level exposure.
  stock_beta_matrix <- matrix(NA_real_, nrow = n_stocks, ncol = k,
                              dimnames = list(stock_id, factor_names))

  for (i in seq_len(n_stocks)) {
    stock_beta_matrix[i, ] <- sector_beta_matrix[stock_sector[i], ] +
      rnorm(k, mean = 0, sd = stock_sd)
  }

  # A simple positive-definite equicorrelation matrix for factor returns.
  corr <- matrix(factor_correlation, nrow = k, ncol = k)
  diag(corr) <- 1
  eig_min <- min(eigen(corr, symmetric = TRUE, only.values = TRUE)$values)
  if (eig_min <= 0) stop("factor_correlation produced a non-positive-definite matrix")

  sds <- factor_sd[factor_names]
  cov_mat <- diag(sds, nrow = k) %*% corr %*% diag(sds, nrow = k)
  z <- matrix(rnorm(n_obs * k), nrow = n_obs, ncol = k)
  factor_matrix <- z %*% chol(cov_mat)
  colnames(factor_matrix) <- factor_names

  residual_sd <- runif(
    n_stocks,
    min = residual_sd_range[1],
    max = residual_sd_range[2]
  )
  names(residual_sd) <- stock_id

  panel_list <- vector("list", n_stocks)
  for (i in seq_len(n_stocks)) {
    mu <- drop(factor_matrix %*% stock_beta_matrix[i, ])
    y <- mu + rnorm(n_obs, mean = 0, sd = residual_sd[i])

    df <- data.frame(
      stock = stock_id[i],
      sector = stock_sector[i],
      t = seq_len(n_obs),
      excess_return = y,
      stringsAsFactors = FALSE
    )
    df[factor_names] <- as.data.frame(factor_matrix)
    panel_list[[i]] <- df
  }

  panel <- do.call(rbind, panel_list)
  rownames(panel) <- NULL

  sector_truth <- data.frame(
    sector = rep(sector_id, each = k),
    factor = rep(factor_names, times = n_sectors),
    beta = as.vector(t(sector_beta_matrix)),
    stringsAsFactors = FALSE
  )

  stock_truth <- data.frame(
    stock = rep(stock_id, each = k),
    sector = rep(stock_sector, each = k),
    factor = rep(factor_names, times = n_stocks),
    beta = as.vector(t(stock_beta_matrix)),
    stringsAsFactors = FALSE
  )

  market_truth <- data.frame(
    factor = factor_names,
    beta = as.numeric(market_mean),
    stringsAsFactors = FALSE
  )

  list(
    panel = panel,
    truth = list(
      market = market_truth,
      sector = sector_truth,
      stock = stock_truth,
      residual_sd = data.frame(stock = stock_id, sigma = residual_sd)
    ),
    metadata = list(
      n_stocks = n_stocks,
      n_sectors = n_sectors,
      n_obs = n_obs,
      factor_names = factor_names,
      seed = seed
    )
  )
}
