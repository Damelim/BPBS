source("BPBS_TP.R")
source("BBE_TP.R")
source("BPS_TP.R")

set.seed(2)
n = 4000
x1 = runif(n, 0, 1)
x2 = runif(n, 0, 1)
xmat = cbind(x1, x2)
x_pred = expand.grid(
  seq(0, 1, length = 150),
  seq(0, 1, length = 150)
)

f = function(x1, x2){
  mnx = 1/3
  sdx = 0.1
  1/3 * (
    dnorm(x1, mnx, sdx) + dnorm(x1, 1-mnx, sdx) +
      1/2*dnorm(x1, 0.5, sdx/2)
  ) * sin(pi*x2)
}

y = f(x1, x2) + rnorm(n, 0, 0.5)
knotmax_TP = max(
  0, min(50, floor(nrow(xmat)^(1/ncol(xmat))) - 4)
)


### Proposed method: BPBS
proposed = BPBS_TP(
  xmat, x_pred, y,
  n_mcmc_sample = 1000, nburnin = 500,
  a_sigma = 0, b_sigma = 0,
  c_lambda = 0.315, initlambda = 1,
  tau_shape1 = 1000/nrow(xmat), tau_shape2 = 1,
  inittau = 0.5,
  tau_grid = seq(0.001, 0.999, length = 100),
  nu = 1/2, initknot = 0, invkappa2 = 1e-7,
  Jpropsigma = 2, knotmax = knotmax_TP,
  saveparams = TRUE, plot_fit = FALSE, seed = 1
)


### Benchmark method: BBE
bbe = BBE_TP(
  xmat, x_pred, y,
  n_mcmc_sample = 1000, nburnin = 500,
  a_g = 0.5, b_g = 0.5*nrow(xmat), initg = nrow(xmat),
  nu = 1/2, initknot = 0, Jpropsigma = 2,
  knotmax = knotmax_TP,
  saveparams = TRUE, plot_fit = FALSE, seed = 1
)


### Benchmark method: BPS
bps = BPS_TP(
  xmat, x_pred, y,
  n_mcmc_sample = 1000, nburnin = 500,
  a_sigma = 0, b_sigma = 0, initsigma2 = 1,
  lambda_proposal_sd = c(1, 1),
  a_lambda = 0.01, b_lambda = 0.01, initlambda = 0.01,
  knots = c(20, 20),
  saveparams = TRUE, plot_fit = FALSE, seed = 1
)


### Posterior-mean surfaces
plot_surface = function(fit, title){
  z_mat = matrix(
    fit$y_pred,
    length(unique(fit$x_pred[,1])),
    length(unique(fit$x_pred[,2]))
  )
  filled.contour(
    unique(fit$x_pred[,1]), unique(fit$x_pred[,2]), z_mat,
    levels = seq(min(fit$y_pred), max(fit$y_pred), length = 20),
    nlevels = 20, col = terrain.colors(22)[1:20],
    xlab = "x1", ylab = "x2", main = title
  )
}

plot_surface(proposed, "BPBS posterior mean")
plot_surface(bbe, "BBE posterior mean")
plot_surface(bps, "BPS posterior mean")


### Parameter traces
par(mfrow = c(2, 3))
plot.ts(proposed$Jvec[,1], main = "BPBS: J1")
plot.ts(proposed$Jvec[,2], main = "BPBS: J2")
plot.ts(log(proposed$lambda), main = "BPBS: log(lambda)")
plot.ts(proposed$tau, main = "BPBS: tau")
plot.ts(sqrt(proposed$sigma2), main = "BPBS: sigma")

par(mfrow = c(2, 2))
plot.ts(bbe$Jvec[,1], main = "BBE: J1")
plot.ts(bbe$Jvec[,2], main = "BBE: J2")
plot.ts(log(bbe$lambda), main = "BBE: log(lambda)")
plot.ts(sqrt(bbe$sigma2), main = "BBE: sigma")

par(mfrow = c(1, 3))
plot.ts(log(bps$lambda[,1]), main = "BPS: log(lambda1)")
plot.ts(log(bps$lambda[,2]), main = "BPS: log(lambda2)")
plot.ts(sqrt(bps$sigma2), main = "BPS: sigma")


prediction_summary = data.frame(
  x1 = x_pred[,1],
  x2 = x_pred[,2],
  truth = f(x_pred[,1], x_pred[,2]),
  BPBS = proposed$y_pred,
  BBE = bbe$y_pred,
  BPS = bps$y_pred,
  BPBS_lower = proposed$lower,
  BPBS_upper = proposed$upper,
  BBE_lower = bbe$lower,
  BBE_upper = bbe$upper,
  BPS_lower = bps$lower,
  BPS_upper = bps$upper
)
prediction_summary
