# BPBS

Code for "Synergizing Roughness Penalization and Basis Selection in Bayesian Spline Regression".

## Files

```BPBS.R``` : main code for the proposed method.

```MatMultInv.cpp``` : Rcpp code for fast matrix multiplication and inversion via RcppArmadillo and RcppEigen.

```basismatrices.R``` : R code for creating basis and penalty matrix.

```execution_example.R``` : an execution example with simulated dataset.


## Dependencies 

Dependencies are the following R packages: ```splines```, ```Rcpp```, ```RcppAradillo```, ```MCMCpack```, ```mvtnorm```, ```invgamma```.


## Arguments

```BPBS``` has the following arguments.

- ```x``` : sorted numeric vector of $x \in [0,1]$.

- ```x_pred``` : sorted numeric vector of out of sample $x \in [0,1]$.

- ```y``` : numeric vector of response, corresponding to $x$.

- ```n_mcmc_sample``` : integer, number of MCMC samples.

- ```nburnin``` : integer, number of burnin samples. It has to be strictly smaller than ```n_mcmc_sample```.

- ```a_sigma``` : shape hyperparameter $\geq 0$ of $\sigma^2$. Default value is 0.

- ```b_sigma``` : scale hyperparameter $\geq 0$ of $\sigma^2$. Default value is 0.

- ```c_lambda``` : positive rate parameter for $\lambda$. Bigger $c_\lambda$ imposes a more informative exponential prior for $\lambda$ toward 0.

- ```initlambda``` : positive initial value for $\lambda$. Default value is 1.

- ```inittau``` : between 0 and 1. initial value for $\tau$. Default value is 0.

- ```nu``` : $\nu$, the geometric hyperparameter for $J$. Bigger $J$ leads to a less informative prior.

- ```initknot``` : initial integer value for interior knots used. Default is 30, which means $J = 34$.

- ```knotmax``` : maximal interior knots to be considered. Has to be smaller than or equal to $n - 4$. Default value is 90.

- ```saveparams``` : boolean, whether to save the MCMC samples of $\sigma^2$, $J$, $\lambda$, $\tau$. Default is FALSE.

- ```plot_fit``` : boolean, whether to plot the out of sample predictions. Default is FALSE.

- ```seed``` : integer, random seed for random number generation.


## Returns

```BPBS``` returns a list of the following terms.

- ```x``` : x in the argument.

- ```y``` : y in the argument.

- ```x_pred``` : x_pred in the argument.

- ```y_pred``` : Estimated pointwise posterior mean curve at x_pred. Has same length as x_pred.

- ```upper``` : Upper 2.5% quantile of the 95% pointwise coverage at x_pred. Has same length as x_pred.

- ```lower``` : Lower 2.5% quantile of the 95% pointwise coverage at x_pred. Has same length as x_pred.

- ```J``` : only returned when saveparams is TRUE. MCMC samples of $J$ getting rid of burn-in.

- ```sigma2``` : only returned when saveparams is TRUE. MCMC samples of $\sigma^2$ getting rid of burn-in.

- ```lambda``` : only returned when saveparams is TRUE. MCMC samples of $\lambda$ getting rid of burn-in.

- ```tau``` : only returned when saveparams is TRUE. MCMC samples of $\tau$ getting rid of burn-in.

