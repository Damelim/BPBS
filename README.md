# BPBS

Code for "Penalty-Induced Basis Exploration for Bayesian Splines".

The repository contains the proposed BPBS method and two benchmark methods,
BBE and BPS, for univariate and tensor-product regression. The tensor-product
functions handle both bivariate (`D = 2`) and trivariate (`D = 3`) predictors.

## Files

`BPBS_1D.R` : proposed method for univariate regression.

`BPBS_TP.R` : proposed method for tensor-product regression (`D = 2` or `3`).

`BBE_1D.R` : BBE benchmark for univariate regression.

`BBE_TP.R` : BBE benchmark for tensor-product regression (`D = 2` or `3`).

`BPS_1D.R` : fixed-dimension Bayesian P-spline benchmark for univariate
regression.

`BPS_TP.R` : fixed-dimension Bayesian P-spline benchmark for tensor-product
regression (`D = 2` or `3`).

`MatMultInv.cpp` : Rcpp code for fast matrix multiplication and inversion via
RcppArmadillo and RcppEigen.

`basismatrices.R` : R code for creating basis and penalty matrices.

`execution_example_1D.R` : univariate example comparing BPBS, BBE, and BPS on
the same simulated data.

`execution_example_TP_2D.R` : 2D tensor-product example comparing BPBS, BBE,
and BPS.

`execution_example_TP_3D.R` : 3D tensor-product example comparing BPBS, BBE,
and BPS and plotting an `x1`--`x2` slice at a fixed value of `x3`.

## Dependencies

The code uses the following R packages: `splines`, `Rcpp`, `RcppArmadillo`,
`RcppEigen`, `MCMCpack`, `mvtnorm`, `mgcv`, `mvnfast`, and `GIGrvg`.

```r
install.packages(c("Rcpp", "RcppArmadillo", "RcppEigen", "MCMCpack",
                   "mvtnorm", "mgcv", "mvnfast", "GIGrvg"))
```

All R and C++ files should be placed in the same directory before sourcing a
method file or running an example.

## Methods

| Function | Basis dimension | Scale update | Saved parameter draws |
|---|---|---|---|
| `BPBS_1D` | sampled `J` | GIG update for `lambda` | `J`, `sigma2`, `lambda`, `tau` |
| `BPBS_TP` | sampled componentwise `Jvec` | GIG update for `lambda` | `Jvec`, `sigma2`, `lambda`, `tau` |
| `BBE_1D` | sampled `J` | inverse-gamma update for `g` (returned as `lambda`) | `J`, `sigma2`, `lambda` |
| `BBE_TP` | sampled componentwise `Jvec` | inverse-gamma update for `g` (returned as `lambda`) | `Jvec`, `sigma2`, `lambda` |
| `BPS_1D` | fixed by `knots` | inverse-gamma update | `sigma2`, `lambda`, `log_joint` |
| `BPS_TP` | fixed by `knots` | componentwise MH updates | `sigma2`, componentwise `lambda` |

The GIG update is specific to BPBS. The current BPBS implementation samples
`lambda` directly from its generalized inverse Gaussian full conditional using
`GIGrvg::rgig()`. BBE and BPS retain their own benchmark-specific updates.

## Model-dimension prior

BPBS and BBE use

$\pi(J) \propto \exp\{-\nu\frac{J}{L}\log(\frac{J}{L})\}$,

where `L = 4` for a univariate cubic B-spline model. For a `D`-dimensional
tensor-product model, `J` is the full tensor-product basis dimension and
`L = 4^D`. The default is `nu = 1/2`; `nu` is not exponentiated by `D`.
Larger positive values of `nu` penalize dimensions above the base model more
strongly.

The model-dimension Metropolis--Hastings steps in BPBS and BBE use normalized,
truncated discrete-normal proposals. Their acceptance ratios include the
reverse-to-forward proposal ratio, which is generally not one near the
boundaries of the model space.

BPS uses a fixed basis dimension and therefore has neither a dimension prior
nor a `J`/`Jvec` MCMC record.

## Common inputs and outputs

All 1D functions take:

* `x` : numeric predictor vector in `[0,1]`.
* `x_pred` : out-of-sample predictor vector in `[0,1]`.
* `y` : numeric response vector corresponding to `x`.
* `n_mcmc_sample` : total number of MCMC draws.
* `nburnin` : number of burn-in draws; it must be strictly smaller than
  `n_mcmc_sample`.
* `saveparams` : whether to return post-burn-in parameter draws.
* `plot_fit` : whether to plot the out-of-sample fit.
* `seed` : random seed.

They return a list containing `x`, `y`, `x_pred`, posterior mean `y_pred`, and
the pointwise 97.5% and 2.5% quantiles `upper` and `lower`.

All TP functions replace `x` with an `n` by `D` matrix `xmat`, and require an
`n_pred` by `D` matrix `x_pred`. They return `xmat`, `y`, `x_pred`, `y_pred`,
`upper`, and `lower`. The same TP functions handle `D = 2` and `D = 3`; no
separate trivariate implementation is required.

Automatic `plot_fit` output is provided for `D = 2`. For `D = 3`, predictions
and credible intervals are returned normally, and `execution_example_TP_3D.R`
shows how to plot a two-dimensional slice at a fixed value of the third
predictor.

## BPBS arguments and returns

In addition to the common inputs, `BPBS_1D` and `BPBS_TP` use:

* `a_sigma`, `b_sigma` : nonnegative hyperparameters for `sigma2`.
* `c_lambda` : positive rate of the exponential prior on `lambda`. A larger
  value places more prior mass near zero. The default is `0.315`.
* `initlambda` : positive initial value for `lambda`. The default is 1.
* `tau_shape1`, `tau_shape2` : positive beta-prior shape parameters for `tau`.
* `inittau` : initial value for `tau`, strictly between 0 and 1.
* `tau_grid` : numeric grid in `(0,1)` used to update `tau`.
* `nu` : positive dimension-prior strength. The default is `1/2`.
* `initknot` : initial number of interior knots. The default is 0.
* `invkappa2` : positive prior-precision term for the intercept.
* `Jpropsigma` : standard deviation of the discrete-normal dimension proposal.
* `knotmax` : maximum number of interior knots considered.

For `BPBS_1D`, `saveparams = TRUE` adds `J`, `sigma2`, `lambda`, and `tau` to
the returned list. For `BPBS_TP`, it adds `Jvec`, `sigma2`, `lambda`, and
`tau`, where `Jvec` has one column per predictor dimension.

For tensor-product BPBS, the default upper bound is

```r
max(0, min(50, floor(nrow(xmat)^(1/ncol(xmat))) - 4))
```

which generalizes the previous square-root rule to `D` dimensions.

## BBE arguments and returns

In addition to the common inputs, `BBE_1D` and `BBE_TP` use:

* `a_g`, `b_g` : shape and scale hyperparameters of the inverse-gamma prior on
  `g`. The defaults are `0.5` and `0.5*n`, respectively.
* `initg` : positive initial value for `g`; the default is `n`.
* `nu` : positive dimension-prior strength. The default is `1/2`.
* `initknot` : initial number of interior knots. In TP models this may be a
  scalar or a length-`D` vector.
* `Jpropsigma` : standard deviation of the discrete-normal dimension proposal.
* `knotmax` : maximum number of interior knots considered.

For `BBE_1D`, `saveparams = TRUE` adds `J`, `sigma2`, and `lambda` to the
returned list. Here `lambda` contains the sampled values denoted by `g` in the
BBE update equations. For `BBE_TP`, the corresponding dimension record is the
componentwise matrix `Jvec`.

The default `BBE_TP` value of `knotmax` is 7 for `D = 3`, matching the supplied
trivariate implementation. For `D = 2`, its default is

```r
max(0, min(50, floor(nrow(xmat)^(1/ncol(xmat))) - 4))
```

The trivariate BBE code first attempts the ordinary matrix inverse. Only when
`(knotmax + 4)^3 > n` and that inverse fails does it retry after adding
`inverse_perturbation * I`; the default perturbation is `1e-7`. This guard does
not alter the ordinary bivariate update or successful trivariate inversions.

## BPS arguments and returns

`BPS_1D` uses:

* `sigma2_0`, `nu_0` : hyperparameters entering the inverse-gamma update for
  `sigma2`.
* `g_shape`, `g_scale` : inverse-gamma prior hyperparameters for the smoothing
  scale.
* `initg`, `initsigma2` : positive initial values.
* `knots` : fixed number of interior knots; the default is 20.

For `BPS_1D`, `saveparams = TRUE` adds `sigma2`, `lambda`, and `log_joint`.
The object returned as `lambda` is the smoothing-scale draw denoted by `g` in
the update equations.

`BPS_TP` uses:

* `a_sigma`, `b_sigma` : nonnegative hyperparameters for `sigma2`.
* `initsigma2` : positive initial value for `sigma2`.
* `lambda_proposal_sd` : positive scalar or length-`D` vector of initial
  random-walk proposal standard deviations. These are adapted during burn-in
  using the original BPS acceptance-rate rule.
* `a_lambda`, `b_lambda` : inverse-gamma prior hyperparameters for each
  componentwise smoothing scale.
* `initlambda` : positive scalar or length-`D` vector of initial values.
* `knots` : scalar or length-`D` vector of fixed interior-knot counts. The
  defaults are 20 per coordinate for `D = 2` and 10 per coordinate for
  `D = 3`, matching the supplied bivariate and trivariate implementations.

For `BPS_TP`, `saveparams = TRUE` adds `sigma2` and `lambda`. Here `lambda` is
an `(n_mcmc_sample - nburnin)` by `D` matrix. BPS does not sample `tau` or the
basis dimension, so it returns neither `tau` nor `J`/`Jvec` histories.
