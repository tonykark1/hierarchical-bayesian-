# Hierarchical Bayesian Factor Lab

Research project testing whether hierarchical Bayesian partial pooling improves equity factor-beta estimation relative to standard OLS, ridge regression, and independent Bayesian regressions.

## Research question

**When does sector-level hierarchical shrinkage improve out-of-sample factor-exposure estimation?**

The initial hierarchy is:

`Market -> Sector -> Stock`

For stock `i`, factor `k`, and sector `s(i)`:

```text
beta[i,k] ~ Normal(mu[s(i),k], tau[s(i),k])
mu[s,k]   ~ Normal(mu_market[k], kappa[k])
```

The project is deliberately designed so that hierarchical Bayes is allowed to lose. The objective is to identify when pooling helps, when it does not, and when a bad hierarchy causes harmful shrinkage.

## v0.1 scope

- ~100 liquid equities across 8-11 sectors
- Weekly returns
- Factors: MKT, SMB, HML, RMW, CMA, MOM
- Four estimators:
  1. OLS
  2. Ridge
  3. Independent Bayesian regression
  4. Hierarchical Bayesian regression
- Simulation-based parameter recovery before real-data estimation
- Rolling out-of-sample evaluation
- Data-scarcity experiment
- Prior and posterior predictive checks

Not in v0.1: HMMs, DCC, wavelets, EVT, portfolio optimization, or latent regimes.

## Core hypotheses

1. Hierarchical Bayes should reduce beta-estimation error when stock-level samples are short or noisy.
2. The advantage of partial pooling should decline as stock-level sample size grows.
3. Benefits should vary by factor because sector membership is more informative for some exposures than others.
4. Misspecified groupings should reduce or reverse the benefit of hierarchical shrinkage.

## Evaluation

Primary metrics:

- beta RMSE / MAE
- correlation between estimated and future realized exposures
- estimate stability across rolling windows
- predictive log score where available
- interval coverage in simulation

Bayesian diagnostics:

- R-hat
- bulk and tail ESS
- divergences
- prior predictive checks
- posterior predictive checks
- parameter recovery
- prior sensitivity

## Planned structure

```text
R/
  simulation.R
  model_ols.R
  model_ridge.R
  model_bayes.R
  model_hierarchical.R
  diagnostics.R
  evaluation.R
  plotting.R
analysis/
  01_simulation.qmd
  02_data.qmd
  03_baselines.qmd
  04_hierarchical.qmd
  05_validation.qmd
  06_oos.qmd
tests/testthat/
models/
output/figures/
output/tables/
output/models/
```

## First milestone

Before downloading market data, the simulation pipeline must show that the estimators recover known parameters and quantify the conditions under which hierarchical pooling improves RMSE over OLS and ridge.

### Phase 1 acceptance criteria

- deterministic simulation under a fixed seed
- known market, sector, and stock-level factor betas
- OLS baseline recovers parameters as sample size grows
- recovery metrics computed at stock-factor level
- sample-size grid ready for 26, 52, 104, 156, 260, and 520 weekly observations
- tests cover dimensions, reproducibility, and parameter joins

## Workflow

This repository uses small research phases. New model complexity is added only after the previous layer passes simulation and out-of-sample checks.
