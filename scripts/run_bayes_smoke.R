source("R/simulation.R")
source("R/model_ols.R")
source("R/model_ridge.R")
source("R/model_bayes.R")
source("R/hierarchy_stress.R")
source("R/evaluation.R")

factors <- c("MKT", "HML", "MOM")

sim <- simulate_hierarchy_scenario(
  scenario = "strong",
  n_stocks = 12L,
  n_sectors = 3L,
  n_obs = 52L,
  factor_names = factors,
  seed = 2201L
)

panel <- sim$panel
truth <- sim$truth$stock

ols <- fit_ols_betas(panel, factors)
ols$estimator <- "OLS"

ridge <- fit_ridge_betas(
  panel,
  factor_names = factors,
  lambda_grid = c(1e-3, 1e-2, 1e-1, 1, 10),
  max_splits = 3L
)
ridge$estimator <- "Ridge"

message("Fitting hierarchical Bayesian smoke model...")
hb_fit <- fit_hierarchical_bayes(
  panel,
  factor_names = factors,
  chains = 2L,
  iter = 800L,
  warmup = 400L,
  cores = 2L,
  seed = 2202L,
  backend = "rstan",
  sample_prior = "no",
  adapt_delta = 0.95,
  max_treedepth = 12L
)

hb <- extract_hierarchical_betas(hb_fit, panel, factors)

summarise_estimator <- function(est) {
  joined <- join_estimates_to_truth(est, truth, match_sector = TRUE)
  overall <- summarise_beta_recovery(joined)
  overall <- overall[overall$factor == "ALL", , drop = FALSE]
  overall$estimator <- unique(est$estimator)
  overall[, c("estimator", "n", "bias", "mae", "rmse")]
}

recovery <- rbind(
  summarise_estimator(ols),
  summarise_estimator(ridge),
  summarise_estimator(hb)
)
rownames(recovery) <- NULL

hb_joined <- join_estimates_to_truth(hb, truth, match_sector = TRUE)
coverage <- summarise_interval_recovery(hb_joined)
coverage <- coverage[coverage$factor == "ALL", , drop = FALSE]
diagnostics <- summarise_brms_diagnostics(hb_fit, max_treedepth = 12L)

diagnostics$scenario <- "strong"
diagnostics$n_stocks <- 12L
diagnostics$n_sectors <- 3L
diagnostics$n_obs <- 52L
diagnostics$n_factors <- length(factors)

ridge_rmse <- recovery$rmse[recovery$estimator == "Ridge"]
hb_rmse <- recovery$rmse[recovery$estimator == "Hierarchical Bayes"]
comparison <- data.frame(
  scenario = "strong",
  n_obs = 52L,
  ridge_rmse = ridge_rmse,
  hierarchical_rmse = hb_rmse,
  hb_minus_ridge = hb_rmse - ridge_rmse,
  hb_relative_improvement = (ridge_rmse - hb_rmse) / ridge_rmse,
  stringsAsFactors = FALSE
)

if (!dir.exists("output")) dir.create("output", recursive = TRUE)
utils::write.csv(recovery, "output/bayes_smoke_recovery.csv", row.names = FALSE)
utils::write.csv(coverage, "output/bayes_smoke_coverage.csv", row.names = FALSE)
utils::write.csv(diagnostics, "output/bayes_smoke_diagnostics.csv", row.names = FALSE)
utils::write.csv(comparison, "output/bayes_smoke_vs_ridge.csv", row.names = FALSE)

cat("\nRecovery summary:\n")
print(recovery, row.names = FALSE)
cat("\nHierarchical posterior interval coverage:\n")
print(coverage, row.names = FALSE)
cat("\nSampler diagnostics:\n")
print(diagnostics, row.names = FALSE)
cat("\nHierarchical Bayes versus ridge (positive relative improvement = HB better):\n")
print(comparison, row.names = FALSE)

# Fail loudly on invalid posterior geometry. Performance itself is reported, not
# used as a CI pass/fail condition, because one seeded smoke run is evidence rather
# than a universal theorem about estimator ranking.
if (is.finite(diagnostics$max_rhat) && diagnostics$max_rhat > 1.05) {
  stop("Bayesian smoke fit failed R-hat diagnostic")
}
if (diagnostics$divergences > 0L) {
  stop("Bayesian smoke fit produced divergent transitions")
}
