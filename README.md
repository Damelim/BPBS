# BPBS

Code for "Penalty-Induced Basis Exploration for Bayesian Splines".

## Files

`BPBS_1D.R` : main code for the proposed method (univariate).

`BPBS_TP.R` : main code for the proposed method (multivariate tensor product).

`MatMultInv.cpp` : Rcpp code for fast matrix multiplication and inversion via RcppArmadillo and RcppEigen.

`basismatrices.R` : R code for creating basis and penalty matrices.

`execution_example_1D.R` : an execution example with simulated data for the univariate method.

`execution_example_TP_2D.R` : a 2D tensor-product example using `BPBS_TP`.

`execution_example_TP_3D.R` : a 3D tensor-product example using the same `BPBS_TP` function.

## Dependencies

The code uses the following R packages: `splines`, `Rcpp`, `RcppArmadillo`, `RcppEigen`, `MCMCpack`, `mvtnorm`, `mgcv`, `mvnfast`, and `GIGrvg`.

```r
install.packages(c("Rcpp", "RcppArmadillo", "RcppEigen", "MCMCpack",
                   "mvtnorm", "mgcv", "mvnfast", "GIGrvg"))
```

The current implementation samples `lambda` directly from its generalized inverse Gaussian full conditional using `GIGrvg::rgig()`.

## Model-dimension prior

The current prior is $\pi(J) \propto \exp\left\{-\nu\frac{J}{L}\log\left(\frac{J}{L}\right) \right\}$,

where `L = 4` for the univariate cubic B-spline model. For a `D`-dimensional tensor-product model, `J` is the full tensor-product basis dimension and `L = 4^D`. The default is `nu = 1/2` in both `BPBS_1D` and `BPBS_TP`; `nu` is not exponentiated by `D`. Larger positive values of `nu` penalize dimensions above the base model more strongly.

The model-dimension Metropolis--Hastings step uses a normalized, truncated discrete-normal proposal. The acceptance ratio therefore includes the reverse-to-forward proposal ratio, which is generally not one near the boundaries of the model space.

## Arguments for BPBS_1D

`BPBS_1D` has the following arguments.

* `x` : sorted numeric vector of predictors in `[0,1]`.

* `x_pred` : sorted numeric vector of out-of-sample predictor values in `[0,1]`.

* `y` : numeric response vector corresponding to `x`.

* `n_mcmc_sample` : number of MCMC draws.

* `nburnin` : number of burn-in draws; it must be strictly smaller than `n_mcmc_sample`.

* `a_sigma` : nonnegative shape hyperparameter for `sigma^2`. The default is 0.

* `b_sigma` : nonnegative scale hyperparameter for `sigma^2`. The default is 0.

* `c_lambda` : positive rate of the exponential prior on `lambda`. A larger value places more prior mass near zero. The default is 0.315.

* `initlambda` : positive initial value for `lambda`. The default is 1.

* `tau_shape1`, `tau_shape2` : positive beta-prior shape parameters for `tau`. The defaults are `1000/n` and 1.

* `inittau` : initial value for `tau`, strictly between 0 and 1. The default is 0.5.

* `tau_grid` : numeric grid in `(0,1)` used to update `tau`.

* `nu` : positive penalty-strength hyperparameter in the dimension prior above. The default is `1/2`; larger values give stronger shrinkage toward the base dimension.

* `initknot` : initial number of interior knots. The default is 0, corresponding to `J = 4`.

* `invkappa2` : positive prior-precision term for the intercept. The default is `1e-7`.

* `Jpropsigma` : standard deviation of the discrete-normal model-dimension proposal. The default is 2.

* `knotmax` : maximum number of interior knots considered. It must be chosen so that the basis matrices used by the sampler are numerically well-defined. The default is 90.

* `saveparams` : whether to return post-burn-in draws of `sigma^2`, `J`, `lambda`, and `tau`. The default is `FALSE`.

* `plot_fit` : whether to plot the out-of-sample fit. The default is `FALSE`.

* `seed` : random seed.

## Returns for BPBS_1D

`BPBS_1D` returns a list with `x`, `y`, `x_pred`, the posterior mean `y_pred`, and the pointwise 97.5% and 2.5% quantiles `upper` and `lower`. If `saveparams = TRUE`, it also returns post-burn-in draws in `J`, `sigma2`, `lambda`, and `tau`.

## Inputs and returns for BPBS_TP

Inputs are analogous to the univariate case, except that `x` is replaced by `xmat`, an `n` by `D` predictor matrix, and `x_pred` is an `n_pred` by `D` prediction matrix. The same `BPBS_TP` implementation handles both bivariate (`D = 2`) and trivariate (`D = 3`) models; no separate trivariate implementation is needed.

For tensor-product models, the default upper bound is

```r
max(0, min(50, floor(nrow(xmat)^(1/ncol(xmat))) - 4))
```

which generalizes the previous square-root rule to `D` dimensions. If `saveparams = TRUE`, `BPBS_TP` returns `Jvec`, an `(n_mcmc_sample - nburnin)` by `D` matrix of componentwise basis dimensions, in addition to `sigma2`, `lambda`, and `tau`.

Automatic `plot_fit` output is provided for `D = 2`. For `D = 3`, predictions and credible intervals are returned normally, and `execution_example_TP_3D.R` illustrates how to plot a two-dimensional slice at a fixed value of the third predictor.
