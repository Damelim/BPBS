library(splines)
library(Rcpp)
library(RcppArmadillo)
library(MCMCpack)
library(mvtnorm)

sourceCpp("MatMultInv.cpp")
source("basismatrices.R")


BBE_1D = function(x, x_pred = seq(0.001, 0.999, by = 0.001), y,
                  n_mcmc_sample = 10000, nburnin = 5000,
                  a_g = 0.5, b_g = 0.5 * length(x), initg = length(x),
                  nu = 1/2, initknot = 0, Jpropsigma = 2,
                  knotmax = 90, saveparams = FALSE,
                  plot_fit = FALSE, seed = 1){

  ### 1. Basic definitions
  n = length(x)
  n_pred = length(x_pred)

  if(length(y) != n){
    stop("length(y) must equal length(x).")
  }
  if(n_mcmc_sample <= nburnin){
    stop("n_mcmc_sample must be larger than nburnin.")
  }

  if(knotmax + 4 > min(n, n_pred)){
    knotmax = min(n, n_pred) - 5
  }
  if(knotmax < 1){
    stop("At least two candidate knot models are required.")
  }

  knotlist = 0:knotmax
  n_models = length(knotlist)
  modelidx = initknot + 1

  if(modelidx < 1 || modelidx > n_models){
    stop("initknot must be between 0 and knotmax.")
  }

  ybar = mean(y)
  ycent = y - ybar
  SSTot = sum(ycent^2)

  if(SSTot <= 0){
    stop("y must have positive sample variance.")
  }

  nkeep = n_mcmc_sample - nburnin
  testpredictions = matrix(NA_real_, nrow = n_pred, ncol = nkeep)
  upper = rep(NA_real_, n_pred)
  lower = rep(NA_real_, n_pred)

  if(saveparams == TRUE){
    sigma2history = rep(NA_real_, nkeep)
    lambdahistory = rep(NA_real_, nkeep)
    Jhistory = rep(NA_real_, nkeep)
  }

  ### 2. Store candidate basis matrices
  Btilde_list = vector("list", n_models)
  BtildetBtilde_list = vector("list", n_models)
  Btilde_pred_list = vector("list", n_models)

  for(k in seq_len(n_models)){
    gbm = get_B_matrix_centered(
      x = x, x_pred = x_pred, degree = 3,
      num_interior_knots = knotlist[k]
    )
    Btilde_list[[k]] = gbm$B[,-1, drop = FALSE]
    Btilde_pred_list[[k]] = gbm$B_pred[,-1, drop = FALSE]
    BtildetBtilde_list[[k]] = eigenMapMatMult(
      t(Btilde_list[[k]]), Btilde_list[[k]]
    )
  }

  spline_order = ncol(Btilde_list[[1]]) + 1
  log_J_prior = function(J_total){
    normalized_J = J_total/spline_order
    -nu * normalized_J * log(normalized_J)
  }

  ### 3. Initial state
  g = initg
  Btilde = Btilde_list[[modelidx]]
  BtildetBtilde = BtildetBtilde_list[[modelidx]]
  p = ncol(Btilde) + 1
  Btildety = crossprod(Btilde, y)
  Covmat = armaInv(BtildetBtilde)
  mn = eigenMapMatMult(Covmat, Btildety)
  SSE = sum((ycent - eigenMapMatMult(Btilde, mn))^2)
  one_Rsq = SSE/SSTot
  logevi = (n-p)/2 * log(1+g) -
    (n-1)/2 * log(1 + g*one_Rsq)

  ### 4. Metropolis-within-Gibbs sampling
  set.seed(seed)

  for(iter in seq_len(n_mcmc_sample)){
    ### 4-1. Sample J
    probvec = dnorm(seq_len(n_models), mean = modelidx, sd = Jpropsigma)
    probvec[modelidx] = 0
    probvec = probvec/sum(probvec)
    modelidx_prop = sample(seq_len(n_models), size = 1,
                           replace = TRUE, prob = probvec)
    q_forward = probvec[modelidx_prop]

    probvec_reverse = dnorm(seq_len(n_models), mean = modelidx_prop,
                            sd = Jpropsigma)
    probvec_reverse[modelidx_prop] = 0
    probvec_reverse = probvec_reverse/sum(probvec_reverse)
    q_reverse = probvec_reverse[modelidx]

    Btilde_prop = Btilde_list[[modelidx_prop]]
    BtildetBtilde_prop = BtildetBtilde_list[[modelidx_prop]]
    p_prop = ncol(Btilde_prop) + 1
    Btildety_prop = crossprod(Btilde_prop, y)
    Covmat_prop = armaInv(BtildetBtilde_prop)
    mn_prop = eigenMapMatMult(Covmat_prop, Btildety_prop)

    SSE_prop = sum((ycent - eigenMapMatMult(Btilde_prop, mn_prop))^2)
    one_Rsq_prop = SSE_prop/SSTot
    logevi_prop = (n-p_prop)/2 * log(1+g) -
      (n-1)/2 * log(1 + g*one_Rsq_prop)

    logpriorratio = log_J_prior(p_prop) - log_J_prior(p)
    logacceptance = log(q_reverse) - log(q_forward) +
      logpriorratio + logevi_prop - logevi

    if(log(runif(1)) <= min(0, logacceptance)){
      modelidx = modelidx_prop
      Btilde = Btilde_prop
      BtildetBtilde = BtildetBtilde_prop
      p = p_prop
      Btildety = Btildety_prop
      Covmat = Covmat_prop
      mn = mn_prop
      SSE = SSE_prop
      one_Rsq = one_Rsq_prop
      logevi = logevi_prop
    }

    ### 4-2. Sample sigma^2
    sigma2 = MCMCpack::rinvgamma(
      n = 1, shape = (n-1)/2,
      scale = SSTot * (1 + g*one_Rsq)/(2 + 2*g)
    )

    ### 4-3. Sample the intercept and centered coefficients
    theta1 = rnorm(1, mean = ybar, sd = sqrt(sigma2/n))
    thetastar = mvtnorm::rmvnorm(
      n = 1, mean = g/(g+1) * mn,
      sigma = sigma2 * g/(g+1) * Covmat
    )

    ### 4-4. Sample g
    g = MCMCpack::rinvgamma(
      n = 1, shape = a_g + (p-1)/2,
      scale = b_g +
        (eigenMapMatMult(
          thetastar,
          eigenMapMatMult(BtildetBtilde/(2*sigma2), t(thetastar))
        ))[1]
    )
    logevi = (n-p)/2 * log(1+g) -
      (n-1)/2 * log(1 + g*one_Rsq)

    ### 4-5. Prediction
    if(iter > nburnin){
      keep = iter - nburnin
      testpredictions[,keep] = theta1 + eigenMapMatMult(
        Btilde_pred_list[[modelidx]], t(thetastar)
      )

      if(saveparams == TRUE){
        sigma2history[keep] = sigma2
        lambdahistory[keep] = g
        Jhistory[keep] = p
      }
    }
  }

  ### 5. Posterior summaries
  model_averaging = rowMeans(testpredictions)
  lower = apply(testpredictions, 1, quantile, probs = 0.025)
  upper = apply(testpredictions, 1, quantile, probs = 0.975)

  if(plot_fit == TRUE){
    plot(x, y, cex = 0.5, xlab = "", ylab = "", main = "")
    title("Data, posterior mean (blue), 95% interval (grey)",
          adj = 0.02, line = -1)
    polygon(c(rev(x_pred), x_pred), c(rev(lower), upper),
            col = adjustcolor("grey", alpha.f = 0.5), border = NA)
    points(x, y, cex = 0.5)
    lines(x_pred, model_averaging, col = "blue", lwd = 1.5)
  }

  out = list(
    "x" = x, "y" = y, "x_pred" = x_pred,
    "y_pred" = model_averaging, "upper" = upper, "lower" = lower
  )

  if(saveparams == TRUE){
    out$J = Jhistory
    out$sigma2 = sigma2history
    out$lambda = lambdahistory
  }

  return(out)
}
