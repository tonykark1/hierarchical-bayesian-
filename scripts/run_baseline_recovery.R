# Run the Phase 1 OLS-vs-ridge parameter-recovery experiment and persist results.

source("R/simulation.R")
source("R/model_ols.R")
source("R/model_ridge.R")
source("R/evaluation.R")

sample_sizes <- c(26L, 52L, 104L, 156L, 260L, 520L)

results <- run_baseline_recovery_experiment(
  sample_sizes = sample_sizes,
  n_stocks = 100L,
  n_sectors = 10L,
  seed = 123L,
  lambda_grid = ridge_lambda_grid(n = 21L, min_lambda = 1e-4, max_lambda = 1e2)
)

dir.create("output", showWarnings = FALSE, recursive = TRUE)
write.csv(results, "output/baseline_recovery.csv", row.names = FALSE)

overall <- results[results$factor == "ALL", , drop = FALSE]
overall <- overall[order(overall$n_obs, overall$estimator), , drop = FALSE]
write.csv(overall, "output/baseline_recovery_overall.csv", row.names = FALSE)

cat("\nOverall beta recovery by estimator and sample size:\n")
print(overall, row.names = FALSE)

# Positive values mean ridge has lower RMSE than OLS.
ols <- overall[overall$estimator == "OLS", c("n_obs", "rmse")]
ridge <- overall[overall$estimator == "Ridge", c("n_obs", "rmse")]
comparison <- merge(ols, ridge, by = "n_obs", suffixes = c("_ols", "_ridge"))
comparison$ridge_rmse_improvement <- comparison$rmse_ols - comparison$rmse_ridge
comparison$ridge_relative_improvement <- comparison$ridge_rmse_improvement / comparison$rmse_ols
write.csv(comparison, "output/baseline_ridge_advantage.csv", row.names = FALSE)

cat("\nRidge RMSE advantage (positive = ridge better):\n")
print(comparison, row.names = FALSE)
