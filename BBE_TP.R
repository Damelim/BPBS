library(splines)
library(Rcpp)
library(RcppArmadillo)
library(MCMCpack)
library(mvnfast)
library(mgcv)

sourceCpp("MatMultInv.cpp")
source("basismatrices.R")


BBE_TP = function(xmat, x_pred, y,
                  n_mcmc_sample = 1000, nburnin = 500,
                  a_g = 0.5, b_g = 0.5 * nrow(xmat), initg = nrow(xmat),
                  nu = 1/2, initknot = 0, Jpropsigma = 2,
                  knotmax = if(ncol(xmat) == 3) 7 else
                    max(0, min(50, floor(nrow(xmat)^(1/ncol(xmat))) - 4)),
                  inverse_perturbation = 1e-7,
                  saveparams = FALSE, plot_fit = FALSE, seed = 1){

  ### 1. Basic definitions
  xmat = as.matrix(xmat)
  x_pred = as.matrix(x_pred)
  D = ncol(xmat)
  n = nrow(xmat)
  n_pred = nrow(x_pred)

  if(!(D %in% c(2, 3))){
    stop("BBE_TP supports two- and three-dimensional predictors only.")
  }
  if(ncol(x_pred) != D){
    stop("x_pred must have the same number of columns as xmat.")
  }
  if(length(y) != n){
    stop("length(y) must equal nrow(xmat).")
  }
  if(n_mcmc_sample <= nburnin){
    stop("n_mcmc_sample must be larger than nburnin.")
  }
  if(knotmax < 1){
    stop("At least two candidate knot models are required.")
  }

  knotlist = 0:knotmax
  n_models = length(knotlist)

  if(length(initknot) == 1){
    initknot = rep(initknot, D)
  }
  if(length(initknot) != D || any(initknot < 0) || any(initknot > knotmax)){
    stop("initknot must be a scalar or a length-D vector between 0 and knotmax.")
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
    Jhistory = matrix(NA_real_, nrow = nkeep, ncol = D)
    colnames(Jhistory) = paste0("J", seq_len(D))
  }

  ### 2. Store componentwise basis matrices
  Btilde_list = vector("list", D)
  Btilde_pred_list = vector("list", D)

  for(d in seq_len(D)){
    Btilde_list[[d]] = vector("list", n_models)
    Btilde_pred_list[[d]] = vector("list", n_models)

    for(k in seq_len(n_models)){
      gbm = get_B_matrix(
        x = xmat[,d], x_pred = x_pred[,d], degree = 3,
        num_interior_knots = knotlist[k]
      )
      Btilde_list[[d]][[k]] = gbm$B
      Btilde_pred_list[[d]][[k]] = gbm$B_pred
    }
  }

  spline_order = ncol(Btilde_list[[1]][[1]]) - knotlist[1]
  J_prior_denominator = spline_order^D
  log_J_prior = function(J_total){
    normalized_J = J_total/J_prior_denominator
    -nu * normalized_J * log(normalized_J)
  }

  # This reproduces the guarded inverse used only by the supplied 3D code.
  use_inverse_perturbation =
    D == 3 && (max(knotlist) + spline_order)^3 > n

  armaInv_selected = function(A){
    if(!use_inverse_perturbation){
      return(armaInv(A))
    }

    tryCatch(
      armaInv(A),
      error = function(e){
        armaInv(A + inverse_perturbation * diag(nrow(A)))
      }
    )
  }

  make_training_state = function(modelidx){
    B_list = lapply(seq_len(D), function(d){
      Btilde_list[[d]][[modelidx[d]]]
    })
    Jvec = vapply(B_list, ncol, numeric(1))
    Bfull = tensor.prod.model.matrix(B_list)
    colmeanvec = colMeans(Bfull)
    colmeanmat = matrix(rep(colmeanvec, n), n,
                        length(colmeanvec), byrow = TRUE)
    Bcentered = Bfull - colmeanmat
    Bcentered[,1] = 1
    Bcentered = Bcentered[,-1, drop = FALSE]
    Btildety = crossprod(Bcentered, y)
    BtildetBtilde = crossprod(Bcentered, Bcentered)
    Covmat = armaInv_selected(BtildetBtilde)
    mn = eigenMapMatMult(Covmat, Btildety)
    SSE = sum((ycent - eigenMapMatMult(Bcentered, mn))^2)

    list(
      B_list = B_list, Jvec = Jvec, J_total = prod(Jvec),
      Btilde = Bcentered, colmeanvec = colmeanvec,
      Btildety = Btildety, BtildetBtilde = BtildetBtilde,
      Covmat = Covmat, mn = mn, SSE = SSE, one_Rsq = SSE/SSTot
    )
  }

  ### 3. Initial state
  modelidx = initknot + 1
  g = initg
  state = make_training_state(modelidx)
  logevi = (n-state$J_total+1)/2 * log(1+g) -
    (n-1)/2 * log(1 + g*state$one_Rsq)

  ### 4. Metropolis-within-Gibbs sampling
  set.seed(seed)

  for(iter in seq_len(n_mcmc_sample)){
    ### 4-1. Sample the basis dimension in each coordinate
    for(d in seq_len(D)){
      probvec = dnorm(seq_len(n_models), mean = modelidx[d],
                      sd = Jpropsigma)
      probvec[modelidx[d]] = 0
      probvec = probvec/sum(probvec)
      modelidx_prop_d = sample(seq_len(n_models), size = 1,
                               replace = TRUE, prob = probvec)
      q_forward = probvec[modelidx_prop_d]

      probvec_reverse = dnorm(seq_len(n_models), mean = modelidx_prop_d,
                              sd = Jpropsigma)
      probvec_reverse[modelidx_prop_d] = 0
      probvec_reverse = probvec_reverse/sum(probvec_reverse)
      q_reverse = probvec_reverse[modelidx[d]]

      modelidx_prop = modelidx
      modelidx_prop[d] = modelidx_prop_d
      state_prop = make_training_state(modelidx_prop)

      logevi_prop = (n-state_prop$J_total+1)/2 * log(1+g) -
        (n-1)/2 * log(1 + g*state_prop$one_Rsq)
      logpriorratio = log_J_prior(state_prop$J_total) -
        log_J_prior(state$J_total)
      logacceptance = logpriorratio + logevi_prop - logevi +
        log(q_reverse) - log(q_forward)

      if(log(runif(1)) <= min(0, logacceptance)){
        modelidx = modelidx_prop
        state = state_prop
        logevi = logevi_prop
      }
    }

    ### 4-2. Sample sigma^2
    sigma2 = MCMCpack::rinvgamma(
      n = 1, shape = (n-1)/2,
      scale = SSTot * (1 + g*state$one_Rsq)/(2 + 2*g)
    )

    ### 4-3. Sample centered coefficients
    thetastar = mvnfast::rmvn(
      n = 1, mu = g/(g+1) * state$mn,
      sigma = sigma2 * g/(g+1) * state$Covmat
    )

    ### 4-4. Sample g
    g = MCMCpack::rinvgamma(
      n = 1, shape = a_g + (state$J_total-1)/2,
      scale = b_g +
        (eigenMapMatMult(
          thetastar,
          eigenMapMatMult(
            state$BtildetBtilde/(2*sigma2), t(thetastar)
          )
        ))[1]
    )
    logevi = (n-state$J_total+1)/2 * log(1+g) -
      (n-1)/2 * log(1 + g*state$one_Rsq)

    ### 4-5. Prediction
    if(iter > nburnin){
      keep = iter - nburnin
      alpha = rnorm(1, mean = ybar, sd = sqrt(sigma2/n))
      B_pred_list = lapply(seq_len(D), function(d){
        Btilde_pred_list[[d]][[modelidx[d]]]
      })
      Btilde_pred = tensor.prod.model.matrix(B_pred_list)
      colmeanmat_pred = matrix(
        rep(state$colmeanvec, n_pred), n_pred,
        length(state$colmeanvec), byrow = TRUE
      )
      Btilde_pred = Btilde_pred - colmeanmat_pred
      Btilde_pred[,1] = 1
      Btilde_pred = Btilde_pred[,-1, drop = FALSE]

      testpredictions[,keep] = alpha + eigenMapMatMult(
        Btilde_pred, t(thetastar)
      )

      if(saveparams == TRUE){
        sigma2history[keep] = sigma2
        lambdahistory[keep] = g
        Jhistory[keep,] = state$Jvec
      }
    }
  }

  ### 5. Posterior summaries
  model_averaging = rowMeans(testpredictions)
  lower = apply(testpredictions, 1, quantile, probs = 0.025)
  upper = apply(testpredictions, 1, quantile, probs = 0.975)

  out = list(
    "xmat" = xmat, "y" = y, "x_pred" = x_pred,
    "y_pred" = model_averaging, "upper" = upper, "lower" = lower
  )

  if(saveparams == TRUE){
    out$Jvec = Jhistory
    out$sigma2 = sigma2history
    out$lambda = lambdahistory
  }

  if(plot_fit == TRUE && D == 2){
    z_mat = matrix(model_averaging, length(unique(x_pred[,1])),
                   length(unique(x_pred[,2])))
    filled.contour(
      unique(x_pred[,1]), unique(x_pred[,2]), z_mat,
      levels = seq(min(model_averaging), max(model_averaging), length = 20),
      nlevels = 20, col = terrain.colors(22)[1:20],
      main = "Test Predictions"
    )
  }

  return(out)
}
