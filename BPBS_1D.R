library(splines)
library(Rcpp)
library(RcppArmadillo)
library(MCMCpack)
library(mvtnorm)
library(GIGrvg)

sourceCpp('MatMultInv.cpp')
source('basismatrices.R')

BPBS_1D = function(x, x_pred = seq(0.001, 0.999, by = 0.001), y,
                   n_mcmc_sample = 10000, nburnin = 5000,
                   a_sigma = 0, b_sigma = 0,
                   c_lambda = 0.315, initlambda = 1,
                   tau_shape1 = 1000/n, tau_shape2 = 1, inittau = 0.5,
                   tau_grid = seq(0.0001, 0.9999, length = 200),
                   nu = 1/2, initknot = 0, invkappa2 = 1e-7, Jpropsigma = 2,
                   knotmax = 90, saveparams = F, plot_fit = F, seed = 1){
  ### 1. basic variable definition.
  tauratio_grid = (1-tau_grid)/tau_grid
  lntaugrid = log(tau_grid)
  n = length(x)
  n_pred = length(x_pred)
  smy = sum(y)
  yty_sumysqnkappa2 = crossprod(y,y) - smy^2 / (n + invkappa2)
  exponent_loglik = 2*a_sigma/2 + n
  addedterm_loglik = 2*b_sigma

  testpredictions = matrix(NA, nrow = n_pred, ncol = n_mcmc_sample - nburnin)
  upper = rep(NA, n_pred); lower = rep(NA, n_pred)

  if(saveparams == TRUE){
    sigma2history = rep(NA, n_mcmc_sample - nburnin)
    lambdahistory = rep(NA, n_mcmc_sample - nburnin)
    tauhistory = rep(NA, n_mcmc_sample - nburnin)
    Jhistory = rep(NA, n_mcmc_sample - nburnin)
  }

  ### 2. Store matrices and eigenvalues for computational efficiency.
  Btilde_list = list(); BtildetBtilde_list = list()
  Btilde_pred_list = list(); Ptilde_list = list()
  knotlist = 0:knotmax

  for(k in 1:length(knotlist)){
    gbm = get_B_matrix_centered(x = x, x_pred = x_pred, degree = 3,
                                num_interior_knots = knotlist[k])
    Btilde_list[[k]] = gbm$B[,-1]
    Btilde_pred_list[[k]] = gbm$B_pred[,-1]
    BtildetBtilde_list[[k]] = eigenMapMatMult(t(Btilde_list[[k]]), Btilde_list[[k]])
    J = k+3
    Ptilde_list[[k]] = (get_P_matrix_centered(n_col = J, deriv = 2))[2:J, 2:J]
  }

  # New dimension prior:
  # pi(J) proportional to exp{-nu * (J/l) * log(J/l)}, where l = 4
  # for cubic B-splines. Deriving l from the base model avoids hard-coding 4.
  spline_order = ncol(Btilde_list[[1]]) + 1
  log_J_prior = function(J_total){
    normalized_J = J_total/spline_order
    -nu * normalized_J * log(normalized_J)
  }

  eigenvalue_list = list()
  negative_eigenvalues = rep(0, length(knotlist))

  for(k in 1:length(knotlist)){
    matt = eigenMapMatMult(armaInv(BtildetBtilde_list[[k]]), Ptilde_list[[k]])
    eigenvalue_list[[k]] = n * Re(eigen(matt, only.values = TRUE)$values)
    eigenvalue_list[[k]] = ifelse(abs(eigenvalue_list[[k]]) < 1e-10,
                                  1e-10, eigenvalue_list[[k]])
    negative_eigenvalues[k] = sum(eigenvalue_list[[k]] < 0)
  }

  if(any(which(negative_eigenvalues > 0)) == TRUE){
    n_models = min(which(negative_eigenvalues > 0)) - 1
  }else{
    n_models = length(BtildetBtilde_list)
  }

  gridsampling_firstterm_list = list()
  for(k in 1:n_models){
    gridsampling_firstterm_list[[k]] =
      0.5 * rowSums(log(1 + outer(tauratio_grid, eigenvalue_list[[k]], FUN = "*"))) +
      0.5 * (k+2) * lntaugrid
  }

  ### 3. initial values
  modelidx = initknot + 1
  lambda = initlambda
  tau = inittau
  Btilde = Btilde_list[[modelidx]]
  BtildetBtilde = BtildetBtilde_list[[modelidx]]
  Btildety = crossprod(Btilde, y)
  Ptilde = Ptilde_list[[modelidx]]
  j_1 = ncol(Btilde); I_j_1 = diag(1, j_1)

  ### 4. Blocked Gibbs Sampling
  set.seed(seed)

  for(iter in 1:n_mcmc_sample){
    Covmat = armaInv((1-tau)/lambda * Ptilde +
                     (n + tau/lambda)*BtildetBtilde/n)
    mn = Covmat %*% Btildety
    SSR = yty_sumysqnkappa2 - crossprod(Btildety, mn); SSR = SSR[1,1]
    logevi = 0.5 * determinant(I_j_1 -
      eigenMapMatMult(Covmat, BtildetBtilde), logarithm = TRUE)$modulus[1] -
      exponent_loglik/2 * log(addedterm_loglik + 0.5*SSR)

    ### 4-1) Sample J (: dimension) from pi(J | lambda, tau, y)
    probvec = dnorm(1:n_models, mean = modelidx, sd = Jpropsigma)
    probvec[modelidx] = 0
    probvec = probvec/sum(probvec)
    modelidx_prop = sample(1:n_models, size = 1, replace = T, prob = probvec)
    q_forward = probvec[modelidx_prop]

    probvec_reverse = dnorm(1:n_models, mean = modelidx_prop, sd = Jpropsigma)
    probvec_reverse[modelidx_prop] = 0
    probvec_reverse = probvec_reverse/sum(probvec_reverse)
    q_reverse = probvec_reverse[modelidx]

    Btilde_prop = Btilde_list[[modelidx_prop]]
    BtildetBtilde_prop = BtildetBtilde_list[[modelidx_prop]]
    Btildety_prop = eigenMapMatMult(t(Btilde_prop), y)
    Ptilde_prop = Ptilde_list[[modelidx_prop]]
    j_1_prop = ncol(Btilde_prop); I_j_1_prop = diag(1, j_1_prop)

    Covmat_prop = armaInv((1-tau)/lambda * Ptilde_prop +
                          (n + tau/lambda)*BtildetBtilde_prop/n)
    mn_prop = Covmat_prop %*% Btildety_prop
    SSR_prop = yty_sumysqnkappa2 - crossprod(Btildety_prop, mn_prop)
    SSR_prop = SSR_prop[1,1]
    logevi_prop = 0.5 * determinant(I_j_1_prop -
      eigenMapMatMult(Covmat_prop, BtildetBtilde_prop),
      logarithm = TRUE)$modulus[1] -
      exponent_loglik/2 * log(addedterm_loglik + 0.5*SSR_prop)

    J_current = j_1 + 1
    J_proposed = j_1_prop + 1
    logpriorratio = log_J_prior(J_proposed) - log_J_prior(J_current)
    logacceptance = logpriorratio + (logevi_prop - logevi) +
      log(q_reverse) - log(q_forward)

    if(log(runif(1)) <= min(0, logacceptance)){
      modelidx = modelidx_prop
      Btilde = Btilde_prop
      BtildetBtilde = BtildetBtilde_prop
      Btildety = Btildety_prop
      Ptilde = Ptilde_prop
      j_1 = j_1_prop; I_j_1 = I_j_1_prop
      Covmat = Covmat_prop
      mn = mn_prop
      SSR = SSR_prop
      logevi = logevi_prop
    }

    ### 4-2) Sample sigma^2 and thetatilde
    sigma_2 = MCMCpack::rinvgamma(n = 1, shape = exponent_loglik/2,
                                  scale = (addedterm_loglik + SSR)/2)
    thetatilde = mvtnorm::rmvnorm(n = 1, mean = mn,
                                  sigma = sigma_2 * Covmat)

    ### 4-3) Sample lambda from its GIG full conditional
    thPth = (thetatilde %*% Ptilde %*% t(thetatilde))[1,1]
    thBtBnth = (thetatilde %*% BtildetBtilde %*% t(thetatilde))[1,1]/n
    gig_b = ((1-tau)*thPth + tau*thBtBnth)/sigma_2
    gig_p = 1 - j_1/2
    lambda = GIGrvg::rgig(n = 1, lambda = gig_p,
                          chi = gig_b, psi = 2*c_lambda)

    ### 4-4) Sample tau via grid sampling
    thPth_lambdasig2 = thPth/(lambda*sigma_2)
    thBtBth_lambdasig2 = thBtBnth/(lambda*sigma_2)
    tauprob_grid = gridsampling_firstterm_list[[modelidx]] -
      0.5 * ((1-tau_grid)*thPth_lambdasig2 +
             tau_grid*thBtBth_lambdasig2)
    tau = sample(tau_grid, size = 1,
                 prob = exp(tauprob_grid - max(tauprob_grid)))

    ### 4-5) Out-of-sample point and interval estimation
    if(iter > nburnin){
      theta1 = rnorm(n = 1, mean = smy/(n + invkappa2),
                     sd = sqrt(sigma_2/(n + invkappa2)))
      testpredictions[,iter - nburnin] = theta1 +
        eigenMapMatMult(Btilde_pred_list[[modelidx]], t(thetatilde))

      if(saveparams == TRUE){
        sigma2history[iter - nburnin] = sigma_2
        lambdahistory[iter - nburnin] = lambda
        tauhistory[iter - nburnin] = tau
        Jhistory[iter - nburnin] = modelidx + 3
      }
    }
  }

  ### 5) posterior mean and credible interval
  model_averaging = rowMeans(testpredictions)
  for(j in 1:n_pred){
    jth_testpredictions = testpredictions[j,]
    upper[j] = quantile(jth_testpredictions, probs = 0.975)
    lower[j] = quantile(jth_testpredictions, probs = 0.025)
  }

  if(plot_fit == TRUE){
    plot(x, y, cex = 0.5, font.main = 2, xlab = "", ylab = "", main = "")
    title("Data (dotted), posterior mean (blue), 95% coverage (grey)",
          adj = 0.02, line = -1)
    polygon(c(rev(x_pred), x_pred), c(rev(lower), upper),
            col = adjustcolor("grey", alpha.f = 0.5), border = NA)
    lines(x_pred, model_averaging, col = "blue", lwd = 1.5)
  }

  out = list("x" = x, "y" = y, "x_pred" = x_pred,
             "y_pred" = model_averaging, "upper" = upper, "lower" = lower)

  if(saveparams == TRUE){
    out$J = Jhistory
    out$sigma2 = sigma2history
    out$lambda = lambdahistory
    out$tau = tauhistory
  }

  return(out)
}
