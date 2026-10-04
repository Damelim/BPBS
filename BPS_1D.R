library(splines)
library(Rcpp)
library(RcppArmadillo)
library(MCMCpack)
library(mvnfast)

sourceCpp("MatMultInv.cpp")
source("basismatrices.R")


BPS_1D = function(x, x_pred = seq(0.001, 0.999, by = 0.001), y,
                  n_mcmc_sample = 10000, nburnin = 5000,
                  sigma2_0 = 0, nu_0 = 0,
                  g_shape = 0.01, g_scale = 0.01, initg = 0.01,
                  initsigma2 = 1, knots = 20,
                  saveparams = FALSE, plot_fit = FALSE, seed = 1){

  ### 1. Basic definitions
  n = length(x)
  n_pred = length(x_pred)

  if(length(y) != n){
    stop("length(y) must equal length(x).")
  }
  if(n_mcmc_sample <= nburnin){
    stop("n_mcmc_sample must be larger than nburnin.")
  }

  y_original = y
  meany = mean(y)
  sdy = sd(y)

  if(!is.finite(sdy) || sdy <= 0){
    stop("y must have positive sample standard deviation.")
  }

  y = (y-meany)/sdy
  nkeep = n_mcmc_sample - nburnin
  testpredictions = matrix(NA_real_, nrow = n_pred, ncol = nkeep)

  if(saveparams == TRUE){
    sigma2history = rep(NA_real_, nkeep)
    lambdahistory = rep(NA_real_, nkeep)
    logjointposteriorhistory = rep(NA_real_, nkeep)
  }

  ### 2. Fixed B-spline basis and penalty
  gbm = get_B_matrix(
    x = x, x_pred = x_pred, degree = 3,
    num_interior_knots = knots
  )
  B = gbm$B
  B_pred = gbm$B_pred
  BtB = eigenMapMatMult(t(B), B)
  p = ncol(BtB)
  P = get_P_matrix(n_col = p, deriv = 2)
  Bty = eigenMapMatMult(t(B), y)

  nu0_n = nu_0 + n
  nu0_sigma20 = nu_0*sigma2_0
  gshape_p2 = g_shape + (p-2)/2

  ### 3. Gibbs sampling
  g = initg
  sigma2 = initsigma2
  set.seed(seed)

  for(iter in seq_len(n_mcmc_sample)){
    ### 3-1. Sample coefficients
    Prec = P/g + BtB/sigma2
    Cov = armaInv(Prec)
    mn = eigenMapMatMult(Cov, Bty/sigma2)
    beta = mvnfast::rmvn(n = 1, mu = mn, sigma = Cov)

    ### 3-2. Sample sigma^2
    fitted_values = eigenMapMatMult(B, t(beta))
    resid = y - fitted_values
    SSR = sum(resid^2)
    sigma2 = MCMCpack::rinvgamma(
      n = 1, shape = nu0_n/2,
      scale = (nu0_sigma20 + SSR)/2
    )

    ### 3-3. Sample g
    thPth = (eigenMapMatMult(
      beta, eigenMapMatMult(P, t(beta))
    ))[1]
    g = MCMCpack::rinvgamma(
      n = 1, shape = gshape_p2,
      scale = g_scale + thPth/2
    )

    ### 3-4. Prediction
    if(iter > nburnin){
      keep = iter - nburnin
      testpredictions[,keep] =
        sdy * eigenMapMatMult(B_pred, t(beta)) + meany

      if(saveparams == TRUE){
        sigma2history[keep] = sigma2
        lambdahistory[keep] = g
        logjointposteriorhistory[keep] =
          -n/2*log(sigma2) -
          sum((y-eigenMapMatMult(B, t(beta)))^2)/(2*sigma2) +
          log(MCMCpack::dinvgamma(g, shape = g_shape, scale = g_scale)) -
          (p-2)/2*log(g) - thPth/(2*g)
      }
    }
  }

  ### 4. Posterior summaries
  model_averaging = rowMeans(testpredictions)
  lower = apply(testpredictions, 1, quantile, probs = 0.025)
  upper = apply(testpredictions, 1, quantile, probs = 0.975)

  if(plot_fit == TRUE){
    plot(x, y_original, cex = 0.5, xlab = "", ylab = "", main = "")
    title("Data, posterior mean (blue), 95% interval (grey)",
          adj = 0.02, line = -1)
    polygon(c(rev(x_pred), x_pred), c(rev(lower), upper),
            col = adjustcolor("grey", alpha.f = 0.5), border = NA)
    points(x, y_original, cex = 0.5)
    lines(x_pred, model_averaging, col = "blue", lwd = 1.5)
  }

  out = list(
    "x" = x, "y" = y_original, "x_pred" = x_pred,
    "y_pred" = model_averaging, "upper" = upper, "lower" = lower
  )

  if(saveparams == TRUE){
    out$sigma2 = sigma2history
    out$lambda = lambdahistory
    out$log_joint = logjointposteriorhistory
  }

  return(out)
}
