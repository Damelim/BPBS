source('BPBS_TP.R')

set.seed(2)
n = 4000
x1 = runif(n, 0, 1)
x2 = runif(n, 0, 1)
x3 = runif(n, 0, 1)
xmat = cbind(x1, x2, x3)

prediction_grid = seq(0, 1, length = 40)
x_pred = expand.grid(prediction_grid, prediction_grid, prediction_grid)

f = function(x1, x2, x3){
  sin(2*pi*x1) + cos(2*pi*x2) + 4*(x3 - 0.5)^2
}

y = f(x1, x2, x3) + rnorm(n, 0, 0.5)

proposed = BPBS_TP(xmat, x_pred, y,
                   n_mcmc_sample = 1000, nburnin = 500,
                   a_sigma = 0, b_sigma = 0,
                   c_lambda = 0.315, initlambda = 1,
                   tau_shape1 = 1000/nrow(xmat), tau_shape2 = 1,
                   inittau = 0.5,
                   tau_grid = seq(0.001, 0.999, length = 100),
                   nu = 1/2, initknot = 0, invkappa2 = 1e-7,
                   Jpropsigma = 2,
                   knotmax = max(0, min(50,
                     floor(nrow(xmat)^(1/ncol(xmat)))-4)),
                   saveparams = T, plot_fit = F, seed = 1)

# Plot the posterior mean surface at the prediction-grid value of x3
# closest to 0.5. BPBS_TP itself plots automatically only when D = 2.
x3_slice = prediction_grid[which.min(abs(prediction_grid - 0.5))]
slice_idx = which(x_pred[,3] == x3_slice)
z_slice = matrix(proposed$y_pred[slice_idx],
                 nrow = length(prediction_grid),
                 ncol = length(prediction_grid))

filled.contour(prediction_grid, prediction_grid, z_slice,
               levels = seq(min(z_slice), max(z_slice), length = 20),
               nlevels = 20, col = terrain.colors(22)[1:20],
               xlab = "x1", ylab = "x2",
               main = paste("Posterior mean at x3 =", round(x3_slice, 3)))

par(mfrow = c(2,3))
plot.ts(proposed$Jvec[,1], main = "J1")
plot.ts(proposed$Jvec[,2], main = "J2")
plot.ts(proposed$Jvec[,3], main = "J3")
plot.ts(log(proposed$lambda), main = "log(lambda)")
plot.ts(proposed$tau, main = "tau")
plot.ts(sqrt(proposed$sigma2), main = "sigma")

DF = data.frame(x1 = proposed$x_pred[,1],
                x2 = proposed$x_pred[,2],
                x3 = proposed$x_pred[,3],
                predictions = proposed$y_pred,
                lower = proposed$lower,
                upper = proposed$upper)
DF
