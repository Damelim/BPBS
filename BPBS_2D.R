library(splines)
library(Rcpp)
library(RcppArmadillo)
library(MCMCpack)
library(mvnfast)
library(mgcv)

sourceCpp('MatMultInv.cpp')
source('basismatrices.R')

BPBS_2D = function(xmat, x_pred, y,
                   n_mcmc_sample = 1000, nburnin = 500,
                   a_sigma = 0, b_sigma = 0,
                   c_lambda = 0.315, initlambda = 1,
                   tau_shape1 = 10000/nrow(xmat), tau_shape2 = 1, inittau = 0.5,
                   tau_grid = seq(0.001, 0.999, length = 100),
                   nu = 100^ncol(xmat), initknot = 0, invkappa2 = 1e-7, Jpropsigma = 2,
                   knotmax = min(50,floor(sqrt(nrow(xmat)))-4), saveparams = F, plot_fit = F, seed = 1){
  
  ### 1. basic variable definition.
  tauratio_grid = (1-tau_grid)/tau_grid  # vector of (1-tau)/tau
  lntaugrid = log(tau_grid) # vector of log(tau).
  
  n = nrow(xmat)
  n_pred = nrow(x_pred)
  
  
  smy = sum(y)
  yty_sumysqnkappa2 = crossprod(y,y) - smy^2 / (n + invkappa2)
  exponent_loglik = 2*a_sigma/2 + n
  addedterm_loglik = 2*b_sigma
  
  
  testpredictions = matrix(NA, nrow = n_pred, ncol = n_mcmc_sample - nburnin)
  upper = rep(NA, n_pred) ; lower = rep(NA, n_pred)
  if(saveparams == TRUE){
    sigma2history = rep(NA, n_mcmc_sample - nburnin) ; lambdahistory = rep(NA, n_mcmc_sample - nburnin) ; tauhistory = rep(NA, n_mcmc_sample - nburnin) ; J1history = rep(NA, n_mcmc_sample - nburnin) ; J2history = rep(NA, n_mcmc_sample - nburnin) 
  }
  
  
  
  ### 2. Store matrices and eigenvalues beforehand for computational efficiency
  
  Btilde_list1 = list() ; Btilde_pred_list1 = list() ; Btilde_list2 = list() ; Btilde_pred_list2 = list()
  
  knotlist = 0:knotmax # from 0 interior knots to upper bound of interior knots
  n_models = length(knotlist)
  for(k in 1:length(knotlist)){
    gbm1 = get_B_matrix(x = xmat[,1], x_pred = x_pred[,1], degree = 3, num_interior_knots = knotlist[k])
    Btilde_list1[[k]] = gbm1$B
    Btilde_pred_list1[[k]] = gbm1$B_pred
    
    gbm2 = get_B_matrix(x = xmat[,2], x_pred = x_pred[,2], degree = 3, num_interior_knots = knotlist[k])
    Btilde_list2[[k]] = gbm2$B
    Btilde_pred_list2[[k]] = gbm2$B_pred
  }
  
    
  ### 3. initial values ###
  modelidx1 = initknot + 1  # initknot: initial number of interior knots to use. We consider 0 interior knots, which is the simplest knot specification.
  modelidx2 = initknot + 1
  lambda = initlambda
  tau = inittau
  
  Btilde1 = Btilde_list1[[modelidx1]] ; J1 = ncol(Btilde1)
  Btilde2 = Btilde_list2[[modelidx2]] ; J2 = ncol(Btilde2)
  Btilde = tensor.prod.model.matrix(list(Btilde1,Btilde2))
  colmeanvec = colMeans(Btilde)
  J = ncol(Btilde)
  colmeanmat = matrix(rep(colmeanvec, n), n, J, byrow = T) 
  Btilde = Btilde - colmeanmat
  Btilde = Btilde[,-1] # get rid of first column
  Btildety = crossprod(Btilde, y)
  BtildetBtilde = crossprod(Btilde,Btilde)
    
  Ptildee = tensor.prod.penalties(list(get_P_matrix(n_col = J1, deriv = 2), get_P_matrix(n_col = J2, deriv = 2)))
  Ptilde1 = Ptildee[[1]][2:(J1*J2), 2:(J1*J2)]
  Ptilde2 = Ptildee[[2]][2:(J1*J2), 2:(J1*J2)]
  Ptilde = Ptilde1 + Ptilde2
  I = diag(1, J1*J2 - 1)
    
    
  ### 4. Blocked Gibbs Sampling ###
  
  set.seed(seed)
    
  for(iter in 1:n_mcmc_sample){
      
    # save current values (including initial)
    Covmat = armaInv( (1-tau)/lambda * Ptilde + (n + tau/lambda)*BtildetBtilde/n ) # covariance for sampling betatilde
    mn = Covmat %*% Btildety # eigenMapMatMult(Covmat, Btildety)
    SSR = yty_sumysqnkappa2 - crossprod(Btildety, mn) ; SSR = SSR[1,1]
    logevi = 0.5 * (determinant(I - eigenMapMatMult(Covmat, BtildetBtilde), logarithm = TRUE)$modulus[1]) - (a_sigma + n/2) * log(b_sigma + 0.5*SSR) 

    ### 4-1. Sample J1
    probvec = dnorm(1:n_models, mean = modelidx1, sd = 2) ; probvec[modelidx1] = 0 
    modelidx1_prop = sample(1:n_models, size = 1, replace = T, prob = probvec)
      
    Btilde1_prop = Btilde_list1[[modelidx1_prop]]
    J1_prop = ncol(Btilde1_prop)
      
    Btilde_prop = tensor.prod.model.matrix(list(Btilde1_prop, Btilde2))
    colmeanvec_prop = colMeans(Btilde_prop)
    J_prop = ncol(Btilde_prop)
    colmeanmat_prop = matrix(rep(colmeanvec_prop, n), n, J_prop, byrow = T)
    Btilde_prop = Btilde_prop - colmeanmat_prop
    Btilde_prop = Btilde_prop[,-1]
    Btildety_prop = crossprod(Btilde_prop, y)
    BtildetBtilde_prop = crossprod(Btilde_prop, Btilde_prop)
      
    Ptildee_prop = tensor.prod.penalties(list(get_P_matrix(n_col = J1_prop, deriv = 2), get_P_matrix(n_col = J2, deriv = 2)))
    Ptilde1_prop = Ptildee_prop[[1]][2:(J1_prop*J2), 2:(J1_prop*J2)]
    Ptilde2_prop = Ptildee_prop[[2]][2:(J1_prop*J2), 2:(J1_prop*J2)]
    Ptilde_prop = Ptilde1_prop + Ptilde2_prop
    I_prop = diag(1, J1_prop * J2 - 1)
      

    Covmat_prop = armaInv( (1-tau)/lambda * Ptilde_prop + (n + tau/lambda)*BtildetBtilde_prop/n )
    mn_prop = Covmat_prop %*% Btildety_prop
    SSR_prop = yty_sumysqnkappa2 - eigenMapMatMult(t(Btildety_prop), mn_prop) ; SSR_prop = SSR_prop[1,1]
    logevi_prop = 0.5 * (determinant(I_prop - eigenMapMatMult(Covmat_prop, BtildetBtilde_prop), logarithm = TRUE)$modulus[1]) - (exponent_loglik)/2 * log(addedterm_loglik + 0.5*SSR_prop)
      
    acceptanceprob = exp( (J_prop-3 + nu/exp(1)) * log(nu / (J_prop-3 + nu/exp(1))) - (J-3 + nu/exp(1)) * log(nu / (J-3 + nu/exp(1))) ) * exp(logevi_prop - logevi)

    if(runif(1) <= acceptanceprob){
        
      modelidx1 = modelidx1_prop # most important!
        
      Btilde1 = Btilde1_prop
      J1 = J1_prop
      J = J_prop
        
      Btilde = Btilde_prop
      Btildety = Btildety_prop
      BtildetBtilde = BtildetBtilde_prop
        
      Ptilde = Ptilde_prop
      I = I_prop
        
      Covmat = Covmat_prop
      mn = mn_prop
      SSR = SSR_prop
      logevi = logevi_prop
    } # end what to do when acceptance
      
      

    ### 4-2. Sample J2
      
    probvec = dnorm(1:n_models, mean = modelidx2, sd = 2) ; probvec[modelidx2] = 0 
    modelidx2_prop = sample(1:n_models, size = 1, replace = T, prob = probvec)
      
    Btilde2_prop = Btilde_list2[[modelidx2_prop]]
    J2_prop = ncol(Btilde2_prop)
      
    Btilde_prop = tensor.prod.model.matrix(list(Btilde1, Btilde2_prop)) 
    colmeanvec_prop = colMeans(Btilde_prop)
    J_prop = ncol(Btilde_prop)
    colmeanmat_prop = matrix(rep(colmeanvec_prop, n), n, J_prop, byrow = T)
    Btilde_prop = Btilde_prop - colmeanmat_prop
    Btilde_prop = Btilde_prop[,-1]
    Btildety_prop = crossprod(Btilde_prop, y)
    BtildetBtilde_prop = crossprod(Btilde_prop, Btilde_prop)
      
    Ptildee_prop = tensor.prod.penalties(list(get_P_matrix(n_col = J1, deriv = 2), get_P_matrix(n_col = J2_prop, deriv = 2)))
    Ptilde1_prop = Ptildee_prop[[1]][2:(J1*J2_prop), 2:(J1*J2_prop)]
    Ptilde2_prop = Ptildee_prop[[2]][2:(J1*J2_prop), 2:(J1*J2_prop)]
    Ptilde_prop = Ptilde1_prop + Ptilde2_prop
    I_prop = diag(1, J1 * J2_prop - 1)
      
    Covmat_prop = armaInv( (1-tau)/lambda * Ptilde_prop + (n + tau/lambda)*BtildetBtilde_prop/n )
    mn_prop = Covmat_prop %*% Btildety_prop
    SSR_prop = yty_sumysqnkappa2 - eigenMapMatMult(t(Btildety_prop), mn_prop) ; SSR_prop = SSR_prop[1,1]
    logevi_prop = 0.5 * (determinant(I_prop - eigenMapMatMult(Covmat_prop, BtildetBtilde_prop), logarithm = TRUE)$modulus[1]) - (exponent_loglik)/2 * log(addedterm_loglik + 0.5*SSR_prop)
      
    acceptanceprob = exp( (J_prop-3 + nu/exp(1)) * log(nu / (J_prop-3 + nu/exp(1))) - (J-3 + nu/exp(1)) * log(nu / (J-3 + nu/exp(1))) ) * exp(logevi_prop - logevi)

    if(runif(1) <= acceptanceprob){
        
      modelidx2 = modelidx2_prop # most important!
        
      Btilde2 = Btilde2_prop
      J2 = J2_prop
      J = J_prop
        
      Btilde = Btilde_prop
      Btildety = Btildety_prop
      BtildetBtilde = BtildetBtilde_prop
        
      Ptilde = Ptilde_prop
      I = I_prop
        
      Covmat = Covmat_prop
      mn = mn_prop
      SSR = SSR_prop
      logevi = logevi_prop
    } # end what to do when acceptance
      
      
    ### 4-3. Sample sigma^2 from p(sigma^2 | J1, J2,J3, g, y)
    sigma_2 = MCMCpack::rinvgamma(n = 1, shape = a_sigma + n/2 , scale = b_sigma + SSR/2) # NOTE THAT sigma_2 IS DIFFERENT from sigma2
      
    ### 4-4. Sample theta1 from p(th1 | J, sigma^2, y) and thetastar from p(thstar | J, sigma^2, g, y)
    thetatilde = mvnfast::rmvn(n = 1, mu = mn, sigma = sigma_2 * Covmat) #, checkSymmetry = F)  
 
    ### 4-5. Sample lambda from p(g | J1,J2, tau, theta, sigma^2, y) USING AUXILIARY VARIABLE GIBBS SAMPLING!
    ### Sample auxiliary variable "hh"
    hh = runif(n = 1, min = 0, max = dexp(x = lambda, rate = c_lambda, log = F))
    ### Sample lambda
    newshape = (J1*J2-1)/2 - 1
    thPth = (thetatilde %*% Ptilde %*% t.default(thetatilde))[1, 1]
    thBtBnth = (thetatilde %*% BtildetBtilde %*% t.default(thetatilde))[1,1] / n
    newscale = ( (1-tau)*thPth + tau*thBtBnth ) /(2*sigma_2)
    UUmax = invgamma::pinvgamma( q = -1/c_lambda * log(hh/c_lambda) , shape = newshape, rate = newscale)
    UU = max(1e-10, runif( n = 1, min = 0, max = UUmax ) )
    lambda = invgamma::qinvgamma(p = UU, shape = newshape, rate = newscale)
      
    ### 4-6. Sample tau from p(tau | J, g, theta, sigma^2, y) using grid sampling
    matt = eigenMapMatMult( armaInv(BtildetBtilde), Ptilde )
    evals = n * Re(eigen(matt, only.values = TRUE)$values)
    evals = ifelse(abs(evals) < 1e-10, 1e-10, evals)
      
    if(sum(evals < 0) > 0){
      #warning("Numeric error: some eigenvalues are negative. Converted to positive")
      evals = abs(evals)
    }
      
    thPth_gsig2 = thPth /(lambda*sigma_2)
    thBtBth_gsig2 = thBtBnth / (lambda*sigma_2)
    tauprob_grid = log(dbeta(tau_grid, tau_shape1, tau_shape2)) + 0.5 * rowSums(log( 1+outer(tauratio_grid, evals, FUN="*" )) ) + 0.5 * (J1*J2-1)*lntaugrid - 0.5 * ((1-tau_grid)*thPth_gsig2 + tau_grid*thBtBth_gsig2)
    tau = sample(tau_grid, size = 1, prob = exp(tauprob_grid - max(tauprob_grid)))
      
    ### 4-7. Out of sample point prediction and interval estimation
    if(iter > nburnin){
      alpha = rnorm(n = 1, mean = smy/(n + invkappa2), sd = sqrt(sigma_2 / (n + invkappa2)))
      Btilde_pred = tensor.prod.model.matrix(list(Btilde_pred_list1[[modelidx1]], Btilde_pred_list2[[modelidx2]]))
      colmeanvec = colMeans(  tensor.prod.model.matrix(list(Btilde1,Btilde2)) )
      colmeanmat = matrix(rep(colmeanvec, n_pred), n_pred, length(colmeanvec), byrow = T)
      Btilde_pred = Btilde_pred - colmeanmat
      Btilde_pred = Btilde_pred[,-1]
        
      tp = eigenMapMatMult(Btilde_pred, t.default(thetatilde))
      testpredictions[,iter - nburnin] = alpha + tp
        
      if(saveparams == TRUE){
        lambdahistory[iter - nburnin] = lambda
        sigma2history[iter - nburnin] = sigma_2
        J1history[iter - nburnin] = modelidx1 + 3
        J2history[iter - nburnin] = modelidx2 + 3
        tauhistory[iter - nburnin] = tau
      }
    } # end if statement
  } # end MCMC
    
  ### 5. posterior mean and credible interval of the estimaed function
  model_averaging = rowMeans(testpredictions) # posterior means of estimated out of sample function values.

  for(j in 1:length(model_averaging)){
    jth_testpredictions = testpredictions[j,]
    
    upper[j] = quantile(jth_testpredictions, probs = 0.975) # upper quantile of pointwise interval estimation of jth test point
    lower[j] = quantile(jth_testpredictions, probs = 0.025) # lower quantile of pointwise interval estimation of jth test point
  }
  
  if(plot_fit == TRUE){
    z_mat = matrix(model_averaging, length(unique(x_pred[,1])), length(unique(x_pred[,2])))
    filled.contour(unique(x_pred[,1]), unique(x_pred[,2]), z_mat, 
                   levels = seq(min(model_averaging),max(model_averaging), length = 20), nlevels = 20,
                    col = terrain.colors(22)[1:20],
                   main = "Test Predictions")

  }
    
  
  if(saveparams == TRUE){
    return(list("xmat" = xmat,
                "y" = y,
                "x_pred" = x_pred,
                "y_pred" = model_averaging,
                "upper" = upper,
                "lower" = lower,
                "J1" = J1history,
                "J2" = J2history,
                "sigma2" = sigma2history,
                "lambda" = lambdahistory,
                "tau" = tauhistory))
  }else{
    return(list("xmat" = xmat,
                "y" = y,
                "x_pred" = x_pred,
                "y_pred" = model_averaging,
                "upper" = upper,
                "lower" = lower))
  }
}





