# Stress-test helpers for the economic hierarchy.
#
# A hierarchical estimator should only receive credit when the grouping contains
# real information. These helpers create DGPs with different hierarchy strength
# and deliberately corrupt sector labels while keeping each stock internally
# assigned to exactly one observed sector.

hierarchy_scenarios <- function() {
  data.frame(
    scenario = c("none", "weak", "balanced", "strong"),
    sector_sd = c(0.00, 0.05, 0.20, 0.30),
    stock_sd = c(0.25, 0.25, 0.15, 0.10),
    description = c(
      "Sector labels contain essentially no beta-location information",
      "Small between-sector signal relative to within-sector dispersion",
      "Default DGP with meaningful sector and stock heterogeneity",
      "Strong sector clustering with relatively tight within-sector dispersion"
    ),
    stringsAsFactors = FALSE
  )
}

simulate_hierarchy_scenario <- function(
    scenario = c("none", "weak", "balanced", "strong"),
    n_stocks = 100L,
    n_sectors = 10L,
    n_obs = 156L,
    seed = 123L,
    ...) {
  scenario <- match.arg(scenario)
  config <- hierarchy_scenarios()
  row <- config[config$scenario == scenario, , drop = FALSE]

  sim <- simulate_factor_panel(
    n_stocks = n_stocks,
    n_sectors = n_sectors,
    n_obs = n_obs,
    sector_sd = row$sector_sd,
    stock_sd = row$stock_sd,
    seed = seed,
    ...
  )
  sim$metadata$hierarchy_scenario <- scenario
  sim$metadata$sector_sd <- row$sector_sd
  sim$metadata$stock_sd <- row$stock_sd
  sim
}

corrupt_sector_labels <- function(panel, fraction = 0.25, seed = 123L) {
  if (!is.numeric(fraction) || length(fraction) != 1L ||
      !is.finite(fraction) || fraction < 0 || fraction > 1) {
    stop("fraction must be a finite number in [0, 1]")
  }
  if (!all(c("stock", "sector") %in% names(panel))) {
    stop("panel must contain stock and sector columns")
  }

  out <- panel
  out$stock <- as.character(out$stock)
  out$sector <- as.character(out$sector)
  out$true_sector <- out$sector

  mapping <- unique(out[, c("stock", "sector"), drop = FALSE])
  counts <- table(mapping$stock)
  if (any(counts != 1L)) {
    stop("Each stock must map to exactly one original sector")
  }

  sectors <- sort(unique(mapping$sector))
  if (fraction > 0 && length(sectors) < 2L) {
    stop("At least two sectors are required to corrupt sector labels")
  }

  n_stocks <- nrow(mapping)
  n_change <- as.integer(round(fraction * n_stocks))
  if (fraction > 0 && n_change == 0L) n_change <- 1L

  set.seed(seed)
  selected <- if (n_change > 0L) sample(seq_len(n_stocks), n_change) else integer()
  observed_sector <- mapping$sector

  for (idx in selected) {
    observed_sector[idx] <- sample(setdiff(sectors, mapping$sector[idx]), 1L)
  }

  new_map <- stats::setNames(observed_sector, mapping$stock)
  out$sector <- unname(new_map[out$stock])
  out$sector_corrupted <- out$sector != out$true_sector

  attr(out, "corruption_fraction_requested") <- fraction
  attr(out, "corruption_fraction_realized") <- mean(
    unique(out[, c("stock", "sector_corrupted")])$sector_corrupted
  )
  out
}

sector_mapping_accuracy <- function(panel) {
  required <- c("stock", "sector", "true_sector")
  missing_cols <- setdiff(required, names(panel))
  if (length(missing_cols) > 0L) {
    stop("panel is missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  mapping <- unique(panel[, required, drop = FALSE])
  data.frame(
    n_stocks = nrow(mapping),
    correct = sum(mapping$sector == mapping$true_sector),
    incorrect = sum(mapping$sector != mapping$true_sector),
    accuracy = mean(mapping$sector == mapping$true_sector),
    stringsAsFactors = FALSE
  )
}
