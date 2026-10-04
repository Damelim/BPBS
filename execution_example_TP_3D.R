source("BPBS_TP.R")
source("BBE_TP.R")
source("BPS_TP.R")

set.seed(2)
n = 4000
x1 = runif(n, 0, 1)
x2 = runif(n, 0, 1)
x3 = runif(n, 0, 1)
xmat = cbind(x1, x2, x3)

prediction_grid = seq(0, 1, length = 40)
x_pred = expand.grid(prediction_grid, prediction_grid, prediction_grid)

f = function(x1, x2, x3){
  sin(2*pi*x1) + cos(2*pi*x2) + 4*(x3-0.5)^2
}

y = f(x1, x2, x3) + rnorm(n, 0, 0.5)
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
# knotmax = 7 reproduces the supplied trivariate BBE implementation.
bbe = BBE_TP(
  xmat, x_pred, y,
  n_mcmc_sample = 1000, nburnin = 500,
  a_g = 0.5, b_g = 0.5*nrow(xmat), initg = nrow(xmat),
  nu = 1/2, initknot = 0, Jpropsigma = 2,
  knotmax = 7, inverse_perturbation = 1e-7,
  saveparams = TRUE, plot_fit = FALSE, seed = 1
)


### Benchmark method: BPS
bps = BPS_TP(
  xmat, x_pred, y,
  n_mcmc_sample = 1000, nburnin = 500,
  a_sigma = 0, b_sigma = 0, initsigma2 = 1,
  lambda_proposal_sd = c(1, 1, 1),
  a_lambda = 0.01, b_lambda = 0.01, initlambda = 0.01,
  knots = c(10, 10, 10),
  saveparams = TRUE, plot_fit = FALSE, seed = 1
)


### Plot x1-x2 posterior-mean slices at x3 closest to 0.5
x3_slice = prediction_grid[which.min(abs(prediction_grid-0.5))]
slice_idx = which(x_pred[,3] == x3_slice)

plot_slice = function(fit, title){
  z_slice = matrix(
    fit$y_pred[slice_idx],
    nrow = length(prediction_grid),
    ncol = length(prediction_grid)
  )
  filled.contour(
    prediction_grid, prediction_grid, z_slice,
    levels = seq(min(z_slice), max(z_slice), length = 20),
    nlevels = 20, col = terrain.colors(22)[1:20],
    xlab = "x1", ylab = "x2",
    main = paste(title, "at x3 =", round(x3_slice, 3))
  )
}

plot_slice(proposed, "BPBS posterior mean")
plot_slice(bbe, "BBE posterior mean")
plot_slice(bps, "BPS posterior mean")


### Parameter traces
par(mfrow = c(2, 3))
plot.ts(proposed$Jvec[,1], main = "BPBS: J1")
plot.ts(proposed$Jvec[,2], main = "BPBS: J2")
plot.ts(proposed$Jvec[,3], main = "BPBS: J3")
plot.ts(log(proposed$lambda), main = "BPBS: log(lambda)")
plot.ts(proposed$tau, main = "BPBS: tau")
plot.ts(sqrt(proposed$sigma2), main = "BPBS: sigma")

par(mfrow = c(2, 3))
plot.ts(bbe$Jvec[,1], main = "BBE: J1")
plot.ts(bbe$Jvec[,2], main = "BBE: J2")
plot.ts(bbe$Jvec[,3], main = "BBE: J3")
plot.ts(log(bbe$lambda), main = "BBE: log(lambda)")
plot.ts(sqrt(bbe$sigma2), main = "BBE: sigma")

par(mfrow = c(2, 2))
plot.ts(log(bps$lambda[,1]), main = "BPS: log(lambda1)")
plot.ts(log(bps$lambda[,2]), main = "BPS: log(lambda2)")
plot.ts(log(bps$lambda[,3]), main = "BPS: log(lambda3)")
plot.ts(sqrt(bps$sigma2), main = "BPS: sigma")


prediction_summary = data.frame(
  x1 = x_pred[,1],
  x2 = x_pred[,2],
  x3 = x_pred[,3],
  truth = f(x_pred[,1], x_pred[,2], x_pred[,3]),
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
