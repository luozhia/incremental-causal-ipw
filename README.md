# Incremental Causal Effect for Time to Treatment Initialization

Code accompanying the paper **"Incremental Causal Effect for Time to Treatment
Initialization"** (Ying, Zhao & Xu, *Proceedings of ICLR 2025*).

The paper studies the *incremental causal effect* of intervening on the
**intensity of time to treatment initiation** — i.e. what happens to an outcome
if the hazard of starting treatment is scaled by a factor $\theta$. It provides
identification **without the positivity assumption**, an inverse-probability-
  weighted (IPW) estimator that reuses standard survival-analysis software, and
root-*n* inference via a multiplier bootstrap.

## The estimand and estimator

For a hazard-ratio intervention $\theta$ applied to the treatment-initiation
hazard $\lambda(t \mid L)$, the incremental causal effect is
$\psi(\theta) = E\{Y_{T(\theta)}\}$. Under
consistency and sequential randomization, it is identified (Theorem 1) and
estimated by reweighting the observed outcomes:
  
  ```
psi_hat(theta) = mean( theta^Delta * exp(-(theta - 1) * cumhaz) * Y )
```

where `Delta` is the event indicator, `cumhaz` is the per-subject cumulative
hazard from a fitted Cox model, and `Y` is the outcome. Pointwise confidence
intervals come from a multiplier (weighted) bootstrap.

## Repository structure

```
incremental-causal-effect/
  ├── README.md
├── R/
  │   ├── simulation.R           # Section 4.1: finite-sample simulation study
│   └── data_application.R     # Section 4.2: rheumatoid arthritis (Methotrexate) analysis
└── .gitignore
```

## `R/simulation.R` — simulation study (Section 4.1)

Evaluates the IPW estimator's finite-sample behaviour: bias, empirical standard
error (SEE), bootstrap standard error (SD), and 95% Wald coverage (CP) across a
grid of $\theta$ values and sample sizes $n = 200, 1000, 5000$. Ground-truth
$\psi(\theta)$ is
approximated with a very large Monte Carlo sample. Output is formatted as LaTeX
table rows (Table 1 in the paper).

> **Reproducing the paper's numbers.** For a quick run the script sets
> `rep_num = 1` (Monte Carlo replications) and `B = 2` (bootstrap draws). The
> published Table 1 uses **R = 1000** replications and **B = 200** bootstrap
> draws; set these to reproduce the results. This is computationally heavy.

## `R/data_application.R` — data application (Section 4.2)

Estimates the incremental effect of the time to initiate **Methotrexate (MTX)**
  on the number of tender joints at one year in a rheumatoid arthritis cohort. The
script starts from an already-prepared, one-row-per-patient data frame `df`,
checks the proportional-hazards assumption (Schoenfeld residuals and cumulative
                                            martingale residuals), drops covariates that violate it, applies the IPW
estimator over $\theta \in [0.2, 5]$ with multiplier-bootstrap pointwise 95%
CIs, and plots $\hat{\psi}(\theta)$.

The header of the script documents the exact columns `df` must contain
(`A`, `Delta`, `Y`, and the baseline covariates in `L`).

> **Data availability.** The patient-level data are **not** included in this
> repository for privacy reasons. The construction of `df` from raw per-visit
> records is omitted; supply your own prepared `df` (as described in the script
                                                     > header) to run the analysis.

## Requirements

R with the following packages:
  
```r
install.packages(c("survival", "ggplot2"))
```

`simulation.R` needs only `survival`. `data_application.R` additionally uses `ggplot2`.

## Citation

```bibtex
@inproceedings{ying2025incremental,
  title     = {Incremental Causal Effect for Time to Treatment Initialization},
  author    = {Ying, Andrew and Zhao, Zhichen and Xu, Ronghui},
  booktitle = {International Conference on Learning Representations (ICLR)},
  year      = {2025}
}
```
