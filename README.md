# BMZDP

R implementation of the BMZ-DP approach introduced in "Bayesian semi-parametric inference for clustered recurrent events with zero-inflation and a terminal event" (Tian et al., 2024, *JRSS-C*).

BMZ-DP is a Bayesian semi-parametric model that jointly characterizes a zero-inflated recurrent event process and a terminal event process for clustered data. A point mass mixture of non-homogeneous Poisson processes describes the recurrent intensity, and shared random effects from different sources bridge the non-terminal and terminal event processes. For robustness, nonparametric Dirichlet processes model the residual of the accelerated failure time model for the survival process as well as the cluster-specific frailty distribution. Posterior inference uses a Markov chain Monte Carlo algorithm.

The package exports one function, `bmzdp()`, which runs the sampler and returns the posterior draws. The baseline hazard of the recurrent process can be piecewise constant (used in the paper's STRIDE analysis and simulation setting 1) or Weibull (simulation setting 2).

## Installation

```r
remotes::install_github("xt83/BMZ-DP")
library(BMZDP)
```

The GitHub repository is called `BMZ-DP`, but R package names cannot contain hyphens, so the package itself is called `BMZDP`. You install from `xt83/BMZ-DP` and load with `library(BMZDP)`.

The only dependency is `stats`, which ships with R.

## Model

For subject *i* in cluster *j*:

- **Zero-inflation.** `y_ij ~ Bernoulli(p_ij)`, with `logit(p_ij) = z_ij' eta`. `y_ij = 1` is the point mass: the subject contributes a likelihood factor of one. `y_ij` is sampled only for subjects with no observed recurrent events (`y_index`). For all other subjects it is fixed at 0.
- **Recurrent process** (when `y_ij = 0`). Intensity `gamma_ij * lambda0(t) * exp(zr_ij' betar + mu_j)`.
  - With `baseline = "piecewise"`, `lambda0(t) = lambda_n[q]` on the interval `(time_q[q], time_q[q+1]]`.
  - With `baseline = "weibull"`, `lambda0(t) = kesi * t^(kesi - 1)`.
- **Terminal process.** Weibull accelerated failure time model with subject-specific shape `sigma_ij`: hazard `sigma_ij t^(sigma_ij - 1) gamma_ij^(-zeta sigma_ij) exp(-(mu_intercept + zd_ij' betad + tao mu_j) sigma_ij)`.
- **Shared random effects.**
  - `log(gamma_ij) ~ N(0, epislon_j)` with `epislon_j ~ InvGamma(shape = a0, scale = b0)`. The parameter `zeta` links `gamma_ij` to the terminal process.
  - The cluster effect `mu_j` enters both processes, scaled by `tao` in the terminal process.
- **Dirichlet processes.** Both are truncated stick-breaking DPs with concentration `alpha` and atom prior `N(0, v2)`.
  - The cluster frailties `mu_j` use `num_dp_cluster` atoms (`b_dp_cluster`, labels `s_label`).
  - The terminal-event shapes `sigma_ij` (the AFT residual) use `num_dp_cluster_n` atoms (`b_dp_cluster_n`, labels `s_label_n`).
- **Priors.**
  - `betar ~ N(0, sigma2r)`, `betad ~ N(0, sigma2d)`, `zeta ~ N(0, sigma2zeta)`, `mu_intercept ~ N(0, sigma2mu_intercept)`, `tao ~ N(0, sigma2tao)`, `eta ~ N(0, sigma2eta)`.
  - Piecewise baseline: `lambda_n[q] ~ N(0, sigma2l)`, kept positive by reflected proposals.
  - Weibull baseline: `kesi ~ Gamma(shape = akesi, scale = bkesi)`.

Sampling steps:

- `epislon`, `y`, the DP labels and the stick-breaking weights are updated by direct (Gibbs-type) draws.
- Everything else uses component-wise random-walk Metropolis updates.
- The proposals for `lambda_n` and for the terminal-event DP atoms `b_dp_cluster_n` are reflected at zero.
- Proposals for `kesi` and `gamma_sub` that are not positive are rejected.

The proposal variance for each `gamma_ij` adapts during sampling. It starts at `metrop_var_gamma_sub`. At every iteration it increases by 5 when that subject's running acceptance rate is above 0.5, and decreases by 10 (floored at 1e-8) when the rate is below 0.2.

## Data format

| Argument | Description |
|---|---|
| `event_times` | List of length *n*. Element *i* is the numeric vector of recurrent event times of subject *i* (`numeric(0)` if none). |
| `time` | Last observed time: terminal event time, or censoring time. |
| `delta` | **Censoring** indicator for the terminal event: `1` = censored (no terminal event observed), `0` = terminal event observed at `time`. |
| `cluster` | Integer cluster codes `1, ..., J` (e.g. hospital or practice). *J* is `max(cluster)`. |
| `zr` | *n* x *pr* covariate matrix of the recurrent process. |
| `zd` | *n* x *pd* covariate matrix of the terminal process. |
| `z` | *n* x *pc* covariate matrix of the logistic model for the zero-inflation probability. Add a column of ones for an intercept. |

The recurrent-event indicator `u_ij` of the paper is derived from `event_times` (`u_ij = 1` if the subject has at least one event).

If your data are in the `observed` matrix layout of the original analysis scripts, convert them as follows:

```r
num_rec <- observed[, pd + pr + 6]
event_times <- lapply(seq_len(nrow(observed)), function(l)
  if (num_rec[l] > 0) observed[l, (pd + pr + 7):(pd + pr + num_rec[l] + 6)] else numeric(0))
time <- observed[, 5]
delta <- observed[, 4]
cluster <- observed[, 2]
zd <- observed[, 6:(pd + 5)]
zr <- observed[, (pd + 6):(pd + pr + 5)]
```

## Example

The data-generating process below is a small version of simulation setting 1 of the paper, with a piecewise constant recurrent baseline.

```r
library(BMZDP)
set.seed(2024)

n <- 60
J <- 6
cluster <- rep(seq_len(J), times = n / J)
zr <- sapply(1:3, function(j) as.numeric(scale(rnorm(n, 0, 0.1))))
zd <- cbind(zr[, 2:3], as.numeric(scale(rnorm(n, 0, 0.1))))
z <- cbind(zr, zd[, 3])

betar <- c(0.4, 0.3, 0.2)
betad <- c(0.2, 0.3, 0.4)
zeta <- 0.1
mu_intercept <- 0.15
tao <- -0.5
kappa <- rgamma(n, shape = 7.5, scale = 1)
y0 <- rbinom(n, 1, 1 / (1 + exp(-drop(z %*% rep(1, 4)))))
gamma <- rlnorm(n, meanlog = 0, sdlog = 0.5)
mu <- rnorm(J, mean = seq(-0.4, 0.4, by = 0.2)[sample.int(5, J, replace = TRUE)], sd = 0.1)

scale_d <- exp(mu_intercept + zeta * log(gamma) + drop(zd %*% betad) + tao * mu[cluster])
death <- rweibull(n, shape = kappa, scale = scale_d)
delta <- rbinom(n, 1, 0.5)
time <- ifelse(delta == 0, death, runif(n, 0, death))

time_q <- seq(0, max(time), length.out = 6)
lambda0 <- c(2, 2.3, 2.1, 2.4, 1.9)
event_times <- vector("list", n)
for (i in seq_len(n)) {
  event_times[[i]] <- numeric(0)
  if (y0[i] == 0) {
    rate <- lambda0 * gamma[i] * exp(sum(betar * zr[i, ]) + mu[cluster[i]])
    t_cur <- rexp(1, max(rate))
    while (t_cur <= time[i]) {
      if (runif(1) < rate[findInterval(t_cur, time_q, left.open = TRUE)] / max(rate)) {
        event_times[[i]] <- c(event_times[[i]], t_cur)
      }
      t_cur <- t_cur + rexp(1, max(rate))
    }
  }
}

fit <- bmzdp(
  event_times = event_times, time = time, delta = delta, cluster = cluster,
  zr = zr, zd = zd, z = z,
  baseline = "piecewise", time_q = time_q,
  mcmc_samples = 30, num_dp_cluster = 10, num_dp_cluster_n = 5, alpha = 1, v2 = 1,
  sigma2r = 100, sigma2d = 100, sigma2l = 100, sigma2zeta = 100,
  sigma2mu_intercept = 100, sigma2tao = 100, sigma2eta = 100, a0 = 0.1, b0 = 0.1,
  metrop_var_betar = 0.05, metrop_var_betad = 0.05, metrop_var_lambda_n = 0.01,
  metrop_var_gamma_sub = 50, metrop_var_eta = 0.64, metrop_var_zeta = 0.01,
  metrop_var_mu_intercept = 0.01, metrop_var_tao = 0.01,
  metrop_var_b_dp_cluster = 0.5, metrop_sd_b_dp_cluster_n = 0.01,
  betar_init = 0, betad_init = 0, lambda_n_init = 1, gamma_init = 1, epislon_init = 1,
  eta_init = 0, zeta_init = 0, mu_intercept_init = 0, tao_init = 0, mu_hos_init = 0,
  sigma_init = 1, b_dp_cluster_init = 0, b_dp_cluster_n_init = 3
)

dim(fit$betar)
colMeans(fit$betar[-(1:10), ])
fit$accept$betar / (30 - 1)
```

To fit the Weibull variant, make three changes to the call above:

- Set `baseline = "weibull"`.
- Drop `time_q`, `sigma2l`, `metrop_var_lambda_n` and `lambda_n_init`.
- Add `akesi = 1, bkesi = 1, metrop_var_kesi = 0.01, kesi_init = 1`.

Thirty iterations are only a smoke test. The paper's analyses used 5,000 (simulations) and 10,000 (STRIDE) iterations. The sampler is plain R, and most Metropolis steps re-evaluate the likelihood of all subjects, so run time is roughly proportional to *n* x `mcmc_samples`. Long chains on data sets with thousands of subjects take hours.

For reproducibility, call `set.seed()` before `bmzdp()`. The function makes its first random draws when it initializes `y`, `s_label` and `s_label_n`, in the same order as the original analysis scripts. Given the same data, settings and random number state, it reproduces the chains of those scripts exactly.

## Arguments

Every argument is required and has no default. The arguments specific to one baseline need to be supplied only when that baseline is chosen. Arguments marked "length 1 or *m*" accept a scalar, which is recycled.

The last column gives the values used in the original analysis scripts: the simulations (S1 piecewise, S2 Weibull) and the STRIDE analysis.

### Data

See [Data format](#data-format): `event_times`, `time`, `delta`, `cluster`, `zr`, `zd`, `z`.

### Baseline hazard of the recurrent process

| Argument | Description | Original scripts |
|---|---|---|
| `baseline` | `"piecewise"` or `"weibull"`. | S1, STRIDE: piecewise. S2: Weibull. |
| `time_q` | Piecewise only. Strictly increasing knots `0 = t_0 < ... < t_K`, giving *K* = `length(time_q) - 1` pieces `(t_{q-1}, t_q]`. Every recurrent event time must lie in `(t_0, t_K]`. | S1: `c(0, 1, 2, 3, 4, 5) * t0 / 5` with `t0 = max(time)`. STRIDE: the same, except the last knot is `t0 + 1`. |

### MCMC and Dirichlet process settings

| Argument | Description | Original scripts |
|---|---|---|
| `mcmc_samples` | Number of stored states, including the initial state in row 1. There are `mcmc_samples - 1` sweeps. | 5000 (S1, S2), 10000 (STRIDE) |
| `num_dp_cluster` | Truncation level (number of atoms) of the DP on the cluster frailty `mu_j`. At least 2. | 10 |
| `num_dp_cluster_n` | Truncation level of the DP on the terminal-event shape `sigma_ij`. At least 2. | 5 |
| `alpha` | DP concentration parameter, shared by both DPs. | 1 |
| `v2` | Variance of the `N(0, v2)` prior on the DP atoms, shared by both DPs. | 1 |

### Prior hyperparameters

| Argument | Description | Original scripts |
|---|---|---|
| `sigma2r` | Prior variance of `betar`. | 100 |
| `sigma2d` | Prior variance of `betad`. | 100 |
| `sigma2l` | Piecewise only. Prior variance of each `lambda_n[q]`. | 100 |
| `akesi`, `bkesi` | Weibull only. Shape and scale of the Gamma prior on `kesi`. | 1, 1 |
| `sigma2zeta` | Prior variance of `zeta`. | 100 |
| `sigma2mu_intercept` | Prior variance of `mu_intercept`. | 100 |
| `sigma2tao` | Prior variance of `tao`. | 100 |
| `sigma2eta` | Prior variance of each `eta` coefficient. | 100 |
| `a0`, `b0` | Shape and scale of the inverse-gamma prior on the cluster-level log-frailty variances `epislon_j`. | 0.1, 0.1 |

### Metropolis proposal tuning

`metrop_var_*` arguments are proposal **variances**. `metrop_sd_b_dp_cluster_n` is a proposal **standard deviation**, because the original code hard-coded an sd for this step.

| Argument | Description | Original scripts |
|---|---|---|
| `metrop_var_betar` | Length 1 or *pr*. | S1, S2: 0.05. STRIDE: `c(0.05, 0.05, 0.05, 0.005, 5e-6, 0.1, 0.1)` |
| `metrop_var_betad` | Length 1 or *pd*. | S1, S2: 0.05. STRIDE: as for `betar` |
| `metrop_var_lambda_n` | Piecewise only. Length 1 or *K*. | S1: 0.01. STRIDE: 0.05 |
| `metrop_var_kesi` | Weibull only. | 0.01 |
| `metrop_var_gamma_sub` | Initial value of the adaptive proposal variance of each `gamma_ij`. Length 1 or *n*. | S1, S2: 50. STRIDE: 30 |
| `metrop_var_eta` | Length 1 or *pc*. | S1, S2: 0.64. STRIDE: sd hard-coded to 1e-4, i.e. `metrop_var_eta = 1e-8` |
| `metrop_var_zeta` | Proposal variance for `zeta`. | 0.01 |
| `metrop_var_mu_intercept` | Proposal variance for `mu_intercept`. | 0.01 |
| `metrop_var_tao` | Proposal variance for `tao`. | 0.01 |
| `metrop_var_b_dp_cluster` | For the frailty DP atoms. Length 1 or `num_dp_cluster`. | 0.5 |
| `metrop_sd_b_dp_cluster_n` | Reflected-proposal sd for the shape DP atoms. Length 1 or `num_dp_cluster_n`. | 0.01 |

### Initial values

In the simulations, the regression parameters, `lambda_n`/`kesi`, `gamma`, `epislon`, `mu_hos` and `sigma` started at their true values.

| Argument | Description | Original scripts |
|---|---|---|
| `betar_init` | Length 1 or *pr*. | STRIDE: 0 |
| `betad_init` | Length 1 or *pd*. | STRIDE: 0 |
| `lambda_n_init` | Piecewise only. Length 1 or *K*. | S1: `c(2, 2.3, 2.1, 2.4, 1.9)`. STRIDE: 1 |
| `kesi_init` | Weibull only. | S2: 1.5 |
| `gamma_init` | Subject frailties `gamma_ij`. Length 1 or *n*. | STRIDE: 1 |
| `epislon_init` | Length 1 or *J*. | S1, S2: 0.25. STRIDE: 1 |
| `eta_init` | Length 1 or *pc*. Also sets the initial `p`. | S1, S2: 1. STRIDE: 0 |
| `zeta_init` | Initial `zeta`. | S1, S2: 0.1. STRIDE: 0 |
| `mu_intercept_init` | Initial `mu_intercept`. | S1, S2: 0.15. STRIDE: 0 |
| `tao_init` | Initial `tao`. | S1, S2: -0.5. STRIDE: 0 |
| `mu_hos_init` | Cluster frailties `mu_j`. Length 1 or *J*. | STRIDE: 0 |
| `sigma_init` | Terminal-event shapes `sigma_ij`. Length 1 or *n*. | STRIDE: 1 |
| `b_dp_cluster_init` | Frailty DP atoms. Length 1 or `num_dp_cluster`. | 0 |
| `b_dp_cluster_n_init` | Shape DP atoms. Length 1 or `num_dp_cluster_n`. | 3 |

A few initial states are not arguments:

- The initial stick-breaking weights are uniform (`1 / num_dp_cluster` and `1 / num_dp_cluster_n`).
- The initial DP labels are drawn uniformly at random.
- The initial `y` for subjects in `y_index` is drawn from `Bernoulli(p)`, with `p` computed from `eta_init`.

## Value

`bmzdp()` returns a named list. Chains have one row per stored state, and row 1 holds the initial values. Discard a burn-in before summarizing.

| Element | Dimension | Content |
|---|---|---|
| `baseline` | 1 | The baseline used. |
| `betar` | `mcmc_samples` x *pr* | Recurrent-process regression coefficients. |
| `lambda_n` | `mcmc_samples` x *K* | Piecewise baseline levels (piecewise only). |
| `kesi` | `mcmc_samples` | Weibull shape of the recurrent baseline (Weibull only). |
| `gamma_sub` | `mcmc_samples` x *n* | Subject-level frailties `gamma_ij`. |
| `epislon` | `mcmc_samples` x *J* | Cluster-level variances of `log(gamma_ij)`. |
| `eta` | `mcmc_samples` x *pc* | Zero-inflation logistic coefficients. |
| `p` | `mcmc_samples` x `length(y_index)` | Zero-inflation probabilities of the subjects in `y_index`. |
| `y` | `mcmc_samples` x *n* | Zero-inflation indicators. Columns not in `y_index` are always 0. |
| `y_index` | vector | Subjects without recurrent events, for whom `y` is sampled. |
| `zeta` | `mcmc_samples` | Effect of `log(gamma_ij)` in the terminal process. |
| `mu_intercept` | `mcmc_samples` | Terminal-process intercept. |
| `betad` | `mcmc_samples` x *pd* | Terminal-process regression coefficients. |
| `tao` | `mcmc_samples` | Effect of the cluster frailty `mu_j` in the terminal process. |
| `mu_hos` | `mcmc_samples` x *J* | Cluster frailties `mu_j`. |
| `sigma` | `mcmc_samples` x *n* | Subject-specific Weibull shapes of the terminal process. |
| `s_label` | `mcmc_samples` x *J* (integer) | DP atom label of each cluster, in `1..num_dp_cluster`. |
| `s_label_n` | `mcmc_samples` x *n* (integer) | DP atom label of each subject, in `1..num_dp_cluster_n`. |
| `b_dp_cluster` | `mcmc_samples` x `num_dp_cluster` | Atoms of the frailty DP. |
| `b_dp_cluster_n` | `mcmc_samples` x `num_dp_cluster_n` | Atoms of the shape DP. |
| `w_dp_cluster` | `mcmc_samples` x `num_dp_cluster` | Stick-breaking weights of the frailty DP. |
| `w_dp_cluster_n` | `mcmc_samples` x `num_dp_cluster_n` | Stick-breaking weights of the shape DP. |
| `accept` | list | Metropolis acceptance counts per parameter (`betar`, `lambda_n` or `kesi`, `gamma_sub`, `eta`, `zeta`, `mu_intercept`, `betad`, `tao`, `b_dp_cluster`, `b_dp_cluster_n`). Divide by `mcmc_samples - 1` for rates. |
| `metrop_var_gamma_sub` | *n* | Final adapted proposal variances of `gamma_ij`. |

The original scripts stored the DP labels as one-hot arrays of dimension `mcmc_samples` x *J* x `num_dp_cluster`. Here they are stored as integer labels. To recover the one-hot form for one iteration, use `diag(num_dp_cluster)[fit$s_label[i, ], ]`.

## Citation

Tian, X., Ciarleglio, M., Cai, J., Greene, E. J., Esserman, D., Li, F., and Zhao, Y. (2024). Bayesian semi-parametric inference for clustered recurrent events with zero-inflation and a terminal event. *Journal of the Royal Statistical Society Series C: Applied Statistics*, 73(3), 598–620. https://doi.org/10.1093/jrsssc/qlae003. arXiv:2202.06636.

```bibtex
@article{tian2024bmzdp,
  title   = {Bayesian semi-parametric inference for clustered recurrent events with zero-inflation and a terminal event},
  author  = {Tian, Xinyuan and Ciarleglio, Maria and Cai, Jiachen and Greene, Erich J. and Esserman, Denise and Li, Fan and Zhao, Yize},
  journal = {Journal of the Royal Statistical Society Series C: Applied Statistics},
  volume  = {73},
  number  = {3},
  pages   = {598--620},
  year    = {2024},
  doi     = {10.1093/jrsssc/qlae003}
}
```

## License

Copyright (c) 2026 Xinyuan Tian. All rights reserved.
