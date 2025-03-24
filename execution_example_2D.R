source('BPBS_2D.R')

set.seed(2) 

n = 4000
x1 = runif(n,0,1)  
x2 = runif(n,0,1)  
xmat = cbind(x1,x2)
x_pred = expand.grid(seq(0,1,length = 150), seq(0,1,length = 150))
f = function(x1,x2){
  mnx = 1/3 ; sdx = 0.1
  1/3*(dnorm(x1,mnx,sdx)+dnorm(x1,1-mnx,sdx)+1/2*dnorm(x1,0.5,sdx/2)) * sin(pi*x2) 
}
y = f(x1,x2) + rnorm(n,0,0.5)

proposed = BPBS_2D(xmat, x_pred, y,
                   n_mcmc_sample = 1000, nburnin = 500,
                   a_sigma = 0, b_sigma = 0,
                   c_lambda = 0.315, initlambda = 1,
                   tau_shape1 = 10000/nrow(xmat), tau_shape2 = 1, inittau = 0.5,
                   tau_grid = seq(0.001, 0.999, length = 100),
                   nu = 100^ncol(xmat), initknot = 0, invkappa2 = 1e-7, Jpropsigma = 2,
                   knotmax = min(50,floor(sqrt(nrow(xmat)))-4), saveparams = T, plot_fit = T, seed = 1)


par(mfrow = c(2,3))
plot.ts(proposed$J1)
plot.ts(proposed$J2)
plot.ts(log(proposed$lambda))
plot.ts(proposed$tau)
plot.ts(sqrt(proposed$sigma2))

DF = data.frame(xgrid = proposed$x_pred, predictions = proposed$y_pred, lower = proposed$lower, upper = proposed$upper)
DF
