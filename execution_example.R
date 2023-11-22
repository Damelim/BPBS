source('BPBS.R')

set.seed(1)
n = 500
x = sort(runif(n = n, 0, 1))
f = function(x){1 + sin(2*pi*x)}
n = length(x)
set.seed(1)
y = f(x) + rnorm(n, 0, sd = 0.5)

proposed = BPBS(x, x_pred = seq(0.001, 0.999,by=0.001), y,
                n_mcmc_sample = 10000, nburnin = 5000, 
                a_sigma = 0, b_sigma = 0, 
                c_lambda = 0.315, initlambda = 1,
                inittau = 0.5, tau_grid = seq(0.001,0.999, length = 200),
                nu = 0.9, initknot = 39, invkappa2 = 1e-7,
                knotmax = 90, saveparams = T, plot_fit = T, seed = 1)

par(mfrow = c(2,2))
plot.ts(proposed$J)
plot.ts(log(proposed$lambda))
plot.ts(proposed$tau)
plot.ts(proposed$sigma2)

DF = data.frame(xgrid = proposed$x_pred, predictions = proposed$y_pred, lower = proposed$lower, upper = proposed$upper)
DF