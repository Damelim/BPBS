library(splines)
library(Rcpp)
library(RcppArmadillo)
library(MCMCpack)
library(mvnfast)
library(mgcv)

sourceCpp("MatMultInv.cpp")
source("basismatrices.R")


BPS_TP = function(xmat, x_pred, y,
                  n_mcmc_sample = 1000, nburnin = 500,
                  a_sigma = 0, b_sigma = 0, initsigma2 = 1,
                  lambda_proposal_sd = rep(1, ncol(xmat)),
                  a_lambda = 0.01, b_lambda = 0.01, initlambda = 0.01,
                  knots = rep(if(ncol(xmat) == 2) 20 else 10, ncol(xmat)),
                  saveparams = FALSE, plot_fit = FALSE, seed = 1){

  ### 1. Basic definitions
  xmat = as.matrix(xmat)
  x_pred = as.matrix(x_pred)
  D = ncol(xmat)
  n = nrow(xmat)
  n_pred = nrow(x_pred)

  if(!(D %in% c(2, 3))){
    stop("BPS_TP supports two- and three-dimensional predictors only.")
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
  if(length(knots) == 1){
    knots = rep(knots, D)
  }
  if(length(knots) != D){
    stop("knots must be a scalar or a length-D vector.")
  }
  if(length(lambda_proposal_sd) == 1){
    lambda_proposal_sd = rep(lambda_proposal_sd, D)
  }
  if(length(lambda_proposal_sd) != D || any(lambda_proposal_sd <= 0)){
    stop("lambda_proposal_sd must be positive and have length one or D.")
  }
  if(length(initlambda) == 1){
    initlambda = rep(initlambda, D)
  }
  if(length(initlambda) != D || any(initlambda <= 0)){
    stop("initlambda must be positive and have length one or D.")
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
    lambdahistory = matrix(NA_real_, nrow = nkeep, ncol = D)
    colnames(lambdahistory) = paste0("lambda", seq_len(D))
  }

  ### 2. Fixed tensor-product basis and penalties
  B_list = vector("list", D)
  B_pred_list = vector("list", D)
  Jvec = rep(NA_real_, D)

  for(d in seq_len(D)){
    gbm = get_B_matrix(
      x = xmat[,d], x_pred = x_pred[,d], degree = 3,
      num_interior_knots = knots[d]
    )
    B_list[[d]] = gbm$B
    B_pred_list[[d]] = gbm$B_pred
    Jvec[d] = ncol(B_list[[d]])
  }

  Bfull = tensor.prod.model.matrix(B_list)
  Bfull_pred = tensor.prod.model.matrix(B_pred_list)
  colmeanvec = colMeans(Bfull)
  J_total = ncol(Bfull)

  colmeanmat = matrix(rep(colmeanvec, n), n, J_total, byrow = TRUE)
  Btilde = Bfull - colmeanmat
  Btilde[,1] = 1

  colmeanmat_pred = matrix(rep(colmeanvec, n_pred), n_pred,
                           J_total, byrow = TRUE)
  Btilde_pred = Bfull_pred - colmeanmat_pred
  Btilde_pred[,1] = 1

  Btilde = Btilde[,-1, drop = FALSE]
  Btilde_pred = Btilde_pred[,-1, drop = FALSE]
  BtildetBtilde = crossprod(Btilde, Btilde)
  Btildety = crossprod(Btilde, y)

  P_component = lapply(seq_len(D), function(d){
    get_P_matrix(n_col = Jvec[d], deriv = 2)
  })
  P_full = tensor.prod.penalties(P_component)
  Ptilde = lapply(P_full, function(Pd){
    Pd[2:J_total, 2:J_total, drop = FALSE]
  })

  ### 3. Initial values and adaptation records
  lambda = initlambda
  sigma_2 = initsigma2
  proposal_sd = lambda_proposal_sd
  acceptance = matrix(0, nrow = n_mcmc_sample, ncol = D)

  adapt_proposal_sd = function(current_sd, acceptance_rate){
    if(acceptance_rate < 0.2){
      if(acceptance_rate < 0.02){
        current_sd/3
      }else if(acceptance_rate < 0.05){
        current_sd/2
      }else if(acceptance_rate < 0.1){
        current_sd/1.5
      }else{
        current_sd/1.25
      }
    }else if(acceptance_rate > 0.4){
      if(acceptance_rate > 0.9){
        current_sd*3
      }else if(acceptance_rate > 0.8){
        current_sd*2
      }else if(acceptance_rate > 0.6){
        current_sd*1.5
      }else{
        current_sd*1.25
      }
    }else{
      current_sd
    }
  }

  ### 4. Metropolis-within-Gibbs sampling
  set.seed(seed)

  for(iter in seq_len(n_mcmc_sample)){
    if(iter <= nburnin && iter %% 100 == 0){
      block = (iter-99):iter
      for(d in seq_len(D)){
        proposal_sd[d] = adapt_proposal_sd(
          proposal_sd[d], mean(acceptance[block,d])
        )
      }
    }

    ### 4-1. Sample coefficients
    penalized_precision = BtildetBtilde
    for(d in seq_len(D)){
      penalized_precision = penalized_precision +
        sigma_2/lambda[d] * Ptilde[[d]]
    }
    Cov = sigma_2 * armaInv(penalized_precision)
    mn = eigenMapMatMult(Cov, Btildety)/sigma_2
    thetatilde = mvnfast::rmvn(n = 1, mu = mn, sigma = Cov)

    ### 4-2. Sample sigma^2
    resid = y - eigenMapMatMult(Btilde, t(thetatilde))
    SSR = sum(resid^2)
    sigma_2 = MCMCpack::rinvgamma(
      n = 1, shape = a_sigma + n/2,
      scale = b_sigma + SSR/2
    )

    ### 4-3. Sample each lambda_d
    for(d in seq_len(D)){
      lambda_prop = rnorm(1, mean = lambda[d], sd = proposal_sd[d])

      if(lambda_prop > 0){
        K_lambda = Reduce(`+`, lapply(seq_len(D), function(j){
          Ptilde[[j]]/lambda[j]
        }))
        lambda_prop_vec = lambda
        lambda_prop_vec[d] = lambda_prop
        K_lambda_prop = Reduce(`+`, lapply(seq_len(D), function(j){
          Ptilde[[j]]/lambda_prop_vec[j]
        }))

        acceptanceprob = exp(
          log(MCMCpack::dinvgamma(lambda_prop, a_lambda, b_lambda)) -
            log(MCMCpack::dinvgamma(lambda[d], a_lambda, b_lambda)) +
            0.5 * (
              determinant(K_lambda_prop, logarithm = TRUE)$modulus[1] -
                determinant(K_lambda, logarithm = TRUE)$modulus[1]
            ) -
            0.5 * (
              thetatilde %*% (K_lambda_prop-K_lambda) %*% t(thetatilde)
            )[1,1]
        )

        if(is.na(acceptanceprob)){
          acceptanceprob = 0
        }

        if(runif(1) <= acceptanceprob){
          lambda[d] = lambda_prop
          acceptance[iter,d] = 1
        }
      }
    }

    ### 4-4. Prediction
    if(iter > nburnin){
      keep = iter - nburnin
      tp = eigenMapMatMult(Btilde_pred, t(thetatilde))
      testpredictions[,keep] = sdy*tp + meany

      if(saveparams == TRUE){
        sigma2history[keep] = sigma_2
        lambdahistory[keep,] = lambda
      }
    }
  }

  ### 5. Posterior summaries
  model_averaging = rowMeans(testpredictions)
  lower = apply(testpredictions, 1, quantile, probs = 0.025)
  upper = apply(testpredictions, 1, quantile, probs = 0.975)

  out = list(
    "xmat" = xmat, "y" = y_original, "x_pred" = x_pred,
    "y_pred" = model_averaging, "upper" = upper, "lower" = lower
  )

  if(saveparams == TRUE){
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
