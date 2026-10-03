library(splines)
library(Rcpp)
library(RcppArmadillo)
library(MCMCpack)
library(mvnfast)
library(mgcv)
library(GIGrvg)

sourceCpp('MatMultInv.cpp')
source('basismatrices.R')

BPBS_TP = function(xmat, x_pred, y,
                   n_mcmc_sample = 1000, nburnin = 500,
                   a_sigma = 0, b_sigma = 0,
                   c_lambda = 0.315, initlambda = 1,
                   tau_shape1 = 1000/nrow(xmat), tau_shape2 = 1, inittau = 0.5,
                   tau_grid = seq(0.001, 0.999, length = 100),
                   nu = 1/2, initknot = 0, invkappa2 = 1e-7, Jpropsigma = 2,
                   knotmax = max(0, min(50,
                     floor(nrow(xmat)^(1/ncol(xmat)))-4)),
                   saveparams = F, plot_fit = F, seed = 1){
  ### 1. compute basics
  D = ncol(xmat)
  n = nrow(xmat)
  n_pred = nrow(x_pred)
  tauratio_grid = (1-tau_grid)/tau_grid
  lntaugrid = log(tau_grid)
  smy = sum(y)
  yty_sumysqnkappa2 = crossprod(y,y) - smy^2/(n + invkappa2)
  exponent_loglik = 2*a_sigma/2 + n
  addedterm_loglik = 2*b_sigma

  testpredictions = matrix(NA, nrow = n_pred,
                           ncol = n_mcmc_sample - nburnin)
  upper = rep(NA, n_pred); lower = rep(NA, n_pred)

  if(saveparams == TRUE){
    sigma2history = rep(NA, n_mcmc_sample - nburnin)
    lambdahistory = rep(NA, n_mcmc_sample - nburnin)
    tauhistory = rep(NA, n_mcmc_sample - nburnin)
    Jhistory = matrix(NA, nrow = n_mcmc_sample - nburnin, ncol = D)
  }

  ### 2. Compute basis and prediction matrices for all dimensions and models
  knotlist = 0:knotmax
  n_models = length(knotlist)
  Btilde_list = vector("list", D)
  Btilde_pred_list = vector("list", D)

  for(d in 1:D){
    Btilde_list[[d]] = vector("list", n_models)
    Btilde_pred_list[[d]] = vector("list", n_models)

    for(k in 1:n_models){
      bm = get_B_matrix(x = xmat[,d], x_pred = x_pred[,d], degree = 3,
                        num_interior_knots = knotlist[k])
      Btilde_list[[d]][[k]] = bm$B
      Btilde_pred_list[[d]][[k]] = bm$B_pred
    }
  }

  # New tensor-product dimension prior:
  # pi(J) proportional to exp{-nu * (J/L) * log(J/L)},
  # L = product of the componentwise base dimensions = l^D.
  spline_order = ncol(Btilde_list[[1]][[1]]) - knotlist[1]
  J_prior_denominator = spline_order^D
  log_J_prior = function(J_total){
    normalized_J = J_total/J_prior_denominator
    -nu * normalized_J * log(normalized_J)
  }

  ### 3. Initial values
  modelidx = rep(initknot + 1, D)
  lambda = initlambda
  tau = inittau
  B_list = lapply(1:D, function(d) Btilde_list[[d]][[modelidx[d]]])
  Jvec = sapply(B_list, ncol)
  Btilde = tensor.prod.model.matrix(B_list)
  colmeanvec_train = colMeans(Btilde)
  colmeanmat = matrix(rep(colmeanvec_train, n), n,
                      length(colmeanvec_train), byrow = TRUE)
  Btilde = Btilde - colmeanmat
  Btilde = Btilde[,-1]
  J = ncol(Btilde)
  Btildety = crossprod(Btilde, y)
  BtildetBtilde = crossprod(Btilde, Btilde)
  P_list = lapply(1:D, function(d) get_P_matrix(n_col = Jvec[d], deriv = 2))
  Ptildee = tensor.prod.penalties(P_list)
  Ptilde = Reduce(`+`, lapply(Ptildee,
    function(m) m[2:(J+1), 2:(J+1)]))
  I = diag(1, J)

  set.seed(seed)

  ### 4. Start MCMC
  for(iter in 1:n_mcmc_sample){
    Covmat = armaInv((1-tau)/lambda * Ptilde +
                     (n + tau/lambda)*BtildetBtilde/n)
    mn = Covmat %*% Btildety
    SSR = yty_sumysqnkappa2 - crossprod(Btildety, mn)[1,1]
    logevi = 0.5 * determinant(I -
      eigenMapMatMult(Covmat, BtildetBtilde), logarithm = TRUE)$modulus[1] -
      exponent_loglik/2 * log(addedterm_loglik + 0.5*SSR)

    ### 4-1. Sample J_d for each dimension
    for(d in 1:D){
      probvec = dnorm(1:n_models, mean = modelidx[d], sd = Jpropsigma)
      probvec[modelidx[d]] = 0
      probvec = probvec/sum(probvec)
      modelidx_prop = sample(1:n_models, size = 1, prob = probvec,
                             replace = TRUE)
      q_forward = probvec[modelidx_prop]

      probvec_reverse = dnorm(1:n_models, mean = modelidx_prop,
                              sd = Jpropsigma)
      probvec_reverse[modelidx_prop] = 0
      probvec_reverse = probvec_reverse/sum(probvec_reverse)
      q_reverse = probvec_reverse[modelidx[d]]

      modelidx_new = modelidx
      modelidx_new[d] = modelidx_prop
      B_list_prop = lapply(1:D,
        function(j) Btilde_list[[j]][[modelidx_new[j]]])
      Jvec_prop = sapply(B_list_prop, ncol)
      Btilde_prop = tensor.prod.model.matrix(B_list_prop)
      colmeanvec_prop = colMeans(Btilde_prop)
      colmeanmat_prop = matrix(rep(colmeanvec_prop, n), n,
                               length(colmeanvec_prop), byrow = TRUE)
      Btilde_prop = Btilde_prop - colmeanmat_prop
      Btilde_prop = Btilde_prop[,-1]
      J_prop = ncol(Btilde_prop)
      Btildety_prop = crossprod(Btilde_prop, y)
      BtildetBtilde_prop = crossprod(Btilde_prop, Btilde_prop)
      P_list_prop = lapply(1:D,
        function(j) get_P_matrix(n_col = Jvec_prop[j], deriv = 2))
      Ptildee_prop = tensor.prod.penalties(P_list_prop)
      Ptilde_prop = Reduce(`+`, lapply(Ptildee_prop,
        function(m) m[2:(J_prop+1), 2:(J_prop+1)]))
      I_prop = diag(1, J_prop)

      Covmat_prop = armaInv((1-tau)/lambda * Ptilde_prop +
                            (n + tau/lambda)*BtildetBtilde_prop/n)
      mn_prop = Covmat_prop %*% Btildety_prop
      SSR_prop = yty_sumysqnkappa2 -
        crossprod(Btildety_prop, mn_prop)[1,1]
      logevi_prop = 0.5 * determinant(I_prop -
        eigenMapMatMult(Covmat_prop, BtildetBtilde_prop),
        logarithm = TRUE)$modulus[1] -
        exponent_loglik/2 * log(addedterm_loglik + 0.5*SSR_prop)

      # J and J_prop omit the centered intercept column, so add one here.
      logpriorratio = log_J_prior(J_prop + 1) - log_J_prior(J + 1)
      logacceptance = logpriorratio + (logevi_prop - logevi) +
        log(q_reverse) - log(q_forward)

      if(log(runif(1)) <= min(0, logacceptance)){
        modelidx = modelidx_new
        Btilde = Btilde_prop
        J = J_prop
        Btildety = Btildety_prop
        BtildetBtilde = BtildetBtilde_prop
        Ptilde = Ptilde_prop
        I = I_prop
        Covmat = Covmat_prop
        mn = mn_prop
        SSR = SSR_prop
        logevi = logevi_prop
        Jvec = Jvec_prop
        colmeanvec_train = colmeanvec_prop
      }
    }

    ### 4-2. Sample sigma^2
    sigma_2 = MCMCpack::rinvgamma(1, shape = a_sigma + n/2,
                                  scale = b_sigma + SSR/2)

    ### 4-3. Sample coefficients
    thetatilde = mvnfast::rmvn(1, mu = mn, sigma = sigma_2 * Covmat)

    ### 4-4. Sample lambda from its GIG full conditional
    thPth = (thetatilde %*% Ptilde %*% t.default(thetatilde))[1,1]
    thBtBnth = (thetatilde %*% BtildetBtilde %*%
      t.default(thetatilde))[1,1]/n
    gig_b = ((1-tau)*thPth + tau*thBtBnth)/sigma_2
    gig_p = 1 - J/2
    lambda = GIGrvg::rgig(n = 1, lambda = gig_p,
                          chi = gig_b, psi = 2*c_lambda)

    ### 4-5. Sample tau using grid sampling
    matt = eigenMapMatMult(armaInv(BtildetBtilde), Ptilde)
    evals = n * Re(eigen(matt, only.values = TRUE)$values)
    evals = ifelse(abs(evals) < 1e-10, 1e-10, evals)

    if(sum(evals < 0) > 0){
      evals = abs(evals)
    }

    thPth_lambdasig2 = thPth/(lambda*sigma_2)
    thBtBth_lambdasig2 = thBtBnth/(lambda*sigma_2)
    tauprob_grid = log(dbeta(tau_grid, tau_shape1, tau_shape2)) +
      0.5 * rowSums(log(1 + outer(tauratio_grid, evals, FUN = "*"))) +
      0.5 * J * lntaugrid -
      0.5 * ((1-tau_grid)*thPth_lambdasig2 +
             tau_grid*thBtBth_lambdasig2)
    tau = sample(tau_grid, size = 1,
                 prob = exp(tauprob_grid - max(tauprob_grid)))

    if(iter > nburnin){
      alpha = rnorm(1, mean = smy/(n + invkappa2),
                    sd = sqrt(sigma_2/(n + invkappa2)))
      B_pred_list = lapply(1:D,
        function(d) Btilde_pred_list[[d]][[modelidx[d]]])
      Btilde_pred = tensor.prod.model.matrix(B_pred_list)
      colmeanmat_pred = matrix(rep(colmeanvec_train, n_pred), n_pred,
                               length(colmeanvec_train), byrow = TRUE)
      Btilde_pred = Btilde_pred - colmeanmat_pred
      Btilde_pred = Btilde_pred[,-1]
      tp = eigenMapMatMult(Btilde_pred, t(thetatilde))
      testpredictions[,iter - nburnin] = alpha + tp

      if(saveparams == TRUE){
        lambdahistory[iter - nburnin] = lambda
        sigma2history[iter - nburnin] = sigma_2
        tauhistory[iter - nburnin] = tau
        Jhistory[iter - nburnin,] = Jvec
      }
    }
  }

  model_averaging = rowMeans(testpredictions)
  for(j in 1:length(model_averaging)){
    upper[j] = quantile(testpredictions[j,], probs = 0.975)
    lower[j] = quantile(testpredictions[j,], probs = 0.025)
  }

  out = list("xmat" = xmat, "y" = y, "x_pred" = x_pred,
             "y_pred" = model_averaging, "upper" = upper, "lower" = lower)

  if(saveparams == TRUE){
    out$Jvec = Jhistory
    out$sigma2 = sigma2history
    out$lambda = lambdahistory
    out$tau = tauhistory
  }

  if(plot_fit == TRUE && D == 2){
    z_mat = matrix(model_averaging, length(unique(x_pred[,1])),
                   length(unique(x_pred[,2])))
    filled.contour(unique(x_pred[,1]), unique(x_pred[,2]), z_mat,
                   levels = seq(min(model_averaging), max(model_averaging),
                                length = 20),
                   nlevels = 20, col = terrain.colors(22)[1:20],
                   main = "Test Predictions")
  }

  return(out)
}
