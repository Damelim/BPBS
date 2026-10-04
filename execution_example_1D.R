source("BPBS_1D.R")
source("BBE_1D.R")
source("BPS_1D.R")

set.seed(1)
n = 500
x = sort(runif(n = n, 0, 1))
x_pred = seq(0.001, 0.999, by = 0.001)
f = function(x){1 + sin(2*pi*x)}

set.seed(1)
y = f(x) + rnorm(n, 0, sd = 0.5)


### Proposed method: BPBS
proposed = BPBS_1D(
  x, x_pred, y,
  n_mcmc_sample = 10000, nburnin = 5000,
  a_sigma = 0, b_sigma = 0,
  c_lambda = 0.315, initlambda = 1,
  tau_shape1 = 1000/n, tau_shape2 = 1, inittau = 0.5,
  tau_grid = seq(0.0001, 0.9999, length = 200),
  nu = 1/2, initknot = 0, invkappa2 = 1e-7,
  Jpropsigma = 2, knotmax = 90,
  saveparams = TRUE, plot_fit = FALSE, seed = 1
)


### Benchmark method: BBE
bbe = BBE_1D(
  x, x_pred, y,
  n_mcmc_sample = 10000, nburnin = 5000,
  a_g = 0.5, b_g = 0.5*n, initg = n,
  nu = 1/2, initknot = 0, Jpropsigma = 2,
  knotmax = 90,
  saveparams = TRUE, plot_fit = FALSE, seed = 1
)


### Benchmark method: BPS
bps = BPS_1D(
  x, x_pred, y,
  n_mcmc_sample = 10000, nburnin = 5000,
  sigma2_0 = 0, nu_0 = 0,
  g_shape = 0.01, g_scale = 0.01, initg = 0.01,
  initsigma2 = 1, knots = 20,
  saveparams = TRUE, plot_fit = FALSE, seed = 1
)


### Posterior-mean comparison
plot(x, y, cex = 0.5, xlab = "x", ylab = "y")
lines(x_pred, proposed$y_pred, col = "blue", lwd = 2)
lines(x_pred, bbe$y_pred, col = "red", lwd = 2, lty = 2)
lines(x_pred, bps$y_pred, col = "darkgreen", lwd = 2, lty = 3)
legend(
  "topright", c("BPBS", "BBE", "BPS"),
  col = c("blue", "red", "darkgreen"),
  lty = c(1, 2, 3), lwd = 2, bty = "n"
)


### Parameter traces
par(mfrow = c(2, 2))
plot.ts(proposed$J, main = "BPBS: J")
plot.ts(log(proposed$lambda), main = "BPBS: log(lambda)")
plot.ts(proposed$tau, main = "BPBS: tau")
plot.ts(sqrt(proposed$sigma2), main = "BPBS: sigma")

par(mfrow = c(1, 3))
plot.ts(bbe$J, main = "BBE: J")
plot.ts(log(bbe$lambda), main = "BBE: log(lambda)")
plot.ts(sqrt(bbe$sigma2), main = "BBE: sigma")

par(mfrow = c(1, 2))
plot.ts(log(bps$lambda), main = "BPS: log(lambda)")
plot.ts(sqrt(bps$sigma2), main = "BPS: sigma")


prediction_summary = data.frame(
  x = x_pred,
  truth = f(x_pred),
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
