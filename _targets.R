library(targets)

tar_option_set(
  packages = character(),
  format = "rds"
)

tar_source("R")

list(
  tar_target(
    sample_sizes,
    c(26L, 52L, 104L, 156L, 260L, 520L)
  ),
  tar_target(
    simulation_156,
    simulate_factor_panel(
      n_stocks = 100L,
      n_sectors = 10L,
      n_obs = 156L,
      seed = 123L
    )
  ),
  tar_target(
    ols_156,
    fit_ols_betas(
      simulation_156$panel,
      simulation_156$metadata$factor_names
    )
  ),
  tar_target(
    ols_156_joined,
    join_estimates_to_truth(ols_156, simulation_156$truth$stock)
  ),
  tar_target(
    ols_156_summary,
    summarise_beta_recovery(ols_156_joined)
  ),
  tar_target(
    ols_sample_size_curve,
    run_ols_recovery_experiment(
      sample_sizes = sample_sizes,
      n_stocks = 100L,
      n_sectors = 10L,
      seed = 123L
    )
  )
)
