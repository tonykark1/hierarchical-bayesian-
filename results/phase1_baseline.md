# Phase 1 baseline results

Run: GitHub Actions `Baseline recovery results` run 35581999367  
Commit: `38ab0b2c7510cf476f08f9fd375b288d2eb47c12`

## Question

Before giving hierarchical Bayes credit for shrinkage, how much can an ordinary leakage-safe ridge estimator improve stock-level factor-beta recovery relative to OLS?

The simulation uses 100 stocks, 10 sectors, six factors (MKT, SMB, HML, RMW, CMA, MOM), and known true stock betas. Ridge lambda is selected with expanding-window validation; predictor scaling is estimated inside each training fold only.

## Overall beta recovery

| Weekly observations | OLS RMSE | Ridge RMSE | Ridge RMSE improvement |
|---:|---:|---:|---:|
| 26 | 0.3965 | 0.3121 | 21.28% |
| 52 | 0.2626 | 0.2492 | 5.10% |
| 104 | 0.1824 | 0.1731 | 5.09% |
| 156 | 0.1547 | 0.1497 | 3.21% |
| 260 | 0.1245 | 0.1238 | 0.58% |
| 520 | 0.0790 | 0.0778 | 1.51% |

## Interpretation

The result establishes a non-Bayesian shrinkage hurdle for Phase 2.

- Shrinkage is especially valuable in the data-poor 26-week setting.
- The incremental value of ridge becomes much smaller as the sample grows.
- Therefore hierarchical Bayes should not be judged merely by whether it beats OLS. It must demonstrate value beyond a tuned generic shrinkage estimator, particularly in the small-sample settings where pooling should matter most.

This is a single seeded simulation design, not a general performance claim. Later experiments will repeat across seeds, hierarchy strengths, factors, and deliberately corrupted group labels.

## Reproduction

Run:

```bash
Rscript scripts/run_baseline_recovery.R
```

The workflow also uploads the full factor-level recovery tables as a GitHub Actions artifact.
