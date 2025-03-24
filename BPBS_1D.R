library(splines)
library(Rcpp)
library(RcppArmadillo)
library(MCMCpack)
library(mvtnorm)

sourceCpp('MatMultInv.cpp')
source('basismatrices.R')

BPBS_1D = function(x, x_pred = seq(0.001, 0.999,by=0.001), y,
                   n_mcmc_sample = 10000, nburnin = 5000, 
                   a_sigma = 0, b_sigma = 0,
                   c_lambda = 0.315, initlambda = 1,
                   tau_shape1 = 10000/n, tau_shape2 = 1, inittau = 0.5,
                   tau_grid = seq(0.0001, 0.9999, length = 200),
                   nu = 100, initknot = 0, invkappa2 = 1e-7, Jpropsigma = 2,
                   knotmax = 90, saveparams = F, plot_fit = F, seed = 1){
  
  
  ### 1. basic variable definition.
  tauratio_grid = (1-tau_grid)/tau_grid  # vector of (1-tau)/tau
  lntaugrid = log(tau_grid) # vector of log(tau).
  
  n = length(x)
  n_pred = length(x_pred)
  
  
  smy = sum(y)
  yty_sumysqnkappa2 = crossprod(y,y) - smy^2 / (n + invkappa2)
  exponent_loglik = 2*a_sigma/2 + n
  addedterm_loglik = 2*b_sigma
  
  
  testpredictions = matrix(NA, nrow = n_pred, ncol = n_mcmc_sample - nburnin) # where we save pointwise estimated functions by MCMC iteration.
  upper = rep(NA, n_pred) ; lower = rep(NA, n_pred) # upper 2.5% quantile and lower 2.5% quantile of estimated function, obtained pointwisely at x_test.
  if(saveparams == TRUE){
    sigma2history = rep(NA, n_mcmc_sample - nburnin) ; lambdahistory = rep(NA, n_mcmc_sample - nburnin) ; tauhistory = rep(NA, n_mcmc_sample - nburnin) ; Jhistory = rep(NA, n_mcmc_sample - nburnin) 
  }
  
  
  
  
  ### 2. Store matrices and eigenvalues for computational efficiency.
  
  Btilde_list = list() ; BtildetBtilde_list = list() ; Btilde_pred_list = list() ; Ptilde_list = list()
  
  knotlist = 0:knotmax # from zero interior knots to upper bound of interior knots. We use cubic B-splines, so, the base model has J = 4.
  for(k in 1:length(knotlist)){
    gbm = get_B_matrix_centered(x = x, x_pred = x_pred, degree = 3, num_interior_knots = knotlist[k]) # sp{1, B2(x) - mean(B2(x)), ..., Bj(x) - mean(Bj(x))}
    Btilde_list[[k]] = gbm$B[,-1]
    Btilde_pred_list[[k]] = gbm$B_pred[,-1]
    BtildetBtilde_list[[k]] = eigenMapMatMult(t(Btilde_list[[k]]), Btilde_list[[k]])
    
    J = k+3
    Ptilde_list[[k]] = (get_P_matrix_centered(n_col = J, deriv = 2))[2:J, 2:J] # list of Ptilde = t(Dtilde) %*% Dtilde.
  }
  
  eigenvalue_list = list() # for grid sampling
  negative_eigenvalues = rep(0, length(knotlist))
  
  for(k in 1:length(knotlist)){
    matt = eigenMapMatMult( armaInv(BtildetBtilde_list[[k]]), Ptilde_list[[k]] )
    eigenvalue_list[[k]] = n * Re(eigen(matt, only.values = TRUE)$values)
    eigenvalue_list[[k]] = ifelse(abs(eigenvalue_list[[k]]) < 1e-10, 1e-10, eigenvalue_list[[k]])
    negative_eigenvalues[k] = sum(eigenvalue_list[[k]] < 0)
  }

  ## exception handling: if numerical error causes eigenvalues become negative, get rid of those cases.
  if(any(which(negative_eigenvalues > 0)) == TRUE){
    n_models = min(which(negative_eigenvalues > 0)) - 1
  }else{
    n_models = length(BtildetBtilde_list)
  }
  
  gridsampling_firstterm_list = list()
  for(k in 1:n_models){
    gridsampling_firstterm_list[[k]] = 0.5 * rowSums(log( 1+outer(tauratio_grid, eigenvalue_list[[k]], FUN="*" )) ) + 0.5 * (k+2)*lntaugrid # (k+2) means "j_1", which is J-1.
  }
  
  

  ### 3. initial values
  modelidx = initknot + 1 # initknot: initial number of interior knots to use. We consider 0 interior knots, which is the simplest knot specification.
  lambda = initlambda
  tau = inittau
  Btilde = Btilde_list[[modelidx]]
  BtildetBtilde = BtildetBtilde_list[[modelidx]]
  Btildety = crossprod(Btilde, y)
  Ptilde = Ptilde_list[[modelidx]]
  j_1 = ncol(Btilde)  ; I_j_1 = diag(1, j_1) # j_1 stands for J-1.
  
  
  
  ### 4. Blocked Gibbs Sampling
  
  set.seed(seed) # seed for random sampling
  
  for(iter in 1:n_mcmc_sample){
    
    # save current values (including inital)
    Covmat = armaInv( (1-tau)/lambda * Ptilde + (n + tau/lambda)*BtildetBtilde/n ) # covariance for sampling thetatilde
    mn = Covmat %*% Btildety
    SSR = yty_sumysqnkappa2 - crossprod(Btildety, mn) ; SSR = SSR[1,1]
    logevi = 0.5 * (determinant(I_j_1 - eigenMapMatMult(Covmat, BtildetBtilde), logarithm = TRUE)$modulus[1]) - (exponent_loglik)/2 * log(addedterm_loglik + 0.5*SSR)

    ### 4-1) Sample J (: dimension) from \pi(J | lambda, tau, y)
    probvec = dnorm(1:n_models, mean = modelidx, sd = Jpropsigma) ; probvec[modelidx] = 0 ; probvec = probvec/sum(probvec) 
    modelidx_prop = sample(1:n_models, size = 1, replace = T, prob = probvec)

    Btilde_prop = Btilde_list[[modelidx_prop]]
    BtildetBtilde_prop = BtildetBtilde_list[[modelidx_prop]]
    Btildety_prop = eigenMapMatMult(t(Btilde_prop), y)
    Ptilde_prop = Ptilde_list[[modelidx_prop]]
    j_1_prop = ncol(Btilde_prop) ; I_j_1_prop = diag(1, j_1_prop)
    
    Covmat_prop = armaInv( (1-tau)/lambda * Ptilde_prop + (n + tau/lambda) * BtildetBtilde_prop/n )
    mn_prop = Covmat_prop %*% Btildety_prop
    SSR_prop = yty_sumysqnkappa2 - crossprod(Btildety_prop, mn_prop) ; SSR_prop = SSR_prop[1,1]
    logevi_prop = 0.5 * (determinant(I_j_1_prop - eigenMapMatMult(Covmat_prop, BtildetBtilde_prop), logarithm = TRUE)$modulus[1]) - (exponent_loglik)/2 * log(addedterm_loglik + 0.5*SSR_prop)
    
    acceptanceprob = exp(   (j_1_prop-3 + nu/exp(1)) * log(nu / (j_1_prop-3 + nu/exp(1))) - (j_1-3 + nu/exp(1)) * log(nu / (j_1-3 + nu/exp(1))) ) * exp(logevi_prop - logevi) # proposal ratio * prior ratio * likelihood ratio
    
    if(runif(1) <= acceptanceprob){
      
      modelidx = modelidx_prop # most important!
      
      Btilde = Btilde_prop
      BtildetBtilde = BtildetBtilde_prop
      Btildety = Btildety_prop
      Ptilde = Ptilde_prop
      j_1 = j_1_prop ; I_j_1 = I_j_1_prop # NOTE) I'm tracking p-1, not p itself!
      
      Covmat = Covmat_prop
      mn = mn_prop
      SSR = SSR_prop
      logevi = logevi_prop
    } # end what to do when acceptance
    
    
    ### 4-2) Sample sigma^2 & thetatilde from \pi(sigma^2, thetatilde | J, lambda, tau, y)
    sigma_2 = MCMCpack::rinvgamma(n = 1, shape = exponent_loglik/2 , scale = (addedterm_loglik + SSR)/2) # NOTE THAT sigma_2 IS DIFFERENT from sigma2
    thetatilde = mvtnorm::rmvnorm(n = 1, mean = mn, sigma = sigma_2 * Covmat) #, checkSymmetry = F)  # beta is a row vector
    
    ### 4-3) Sample lambda from p(lambda | J, tau, theta, sigma^2, y) via slice sampling (auxiliary variable MCMC)

    #### 4-3-1) sample auxiliary variable "hh"
    hh = runif(n = 1, min = 0, max = dexp(x = lambda, rate = c_lambda, log = F))
    
    #### 4-3-2) sample lambda
    newshape = j_1/2 - 1   
    thPth = (thetatilde %*% Ptilde %*% t(thetatilde))[1, 1]
    thBtBnth = (thetatilde %*% BtildetBtilde %*% t(thetatilde))[1,1] / n
    newscale = ( (1-tau)*thPth + tau*thBtBnth ) /(2*sigma_2)
    
    UUmax = invgamma::pinvgamma( q = -1/c_lambda * log(hh/c_lambda) , shape = newshape, rate = newscale) # 
    UU = max(1e-10, runif( n = 1, min = 0, max = UUmax ) )
    lambda = invgamma::qinvgamma(p = UU, shape = newshape, rate = newscale)
  
    
    ### 4-4) Sample \tau from \pi(\tau | J,\lambda,\sigma^2, \theta_1, \theta_tilde, y) via grid sampling
    thPth_lambdasig2 = thPth /(lambda*sigma_2)
    thBtBth_lambdasig2 = thBtBnth / (lambda*sigma_2) 
    tauprob_grid = gridsampling_firstterm_list[[modelidx]] - 0.5 * ((1-tau_grid)*thPth_lambdasig2 + tau_grid*thBtBth_lambdasig2) 
    tau = sample(tau_grid, size = 1, prob = exp(tauprob_grid - max(tauprob_grid))) # subtract max(tauprob_grid) before exponentiating for numerical stability.
    
    
    ### 4-5) Out of sample point and interval estimation.
    if(iter > nburnin){
      theta1 = rnorm(n = 1, mean = smy/(n + invkappa2), sd = sqrt(sigma_2 / (n + invkappa2))) # the intercept, theta1.
      testpredictions[,iter - nburnin] = theta1 + eigenMapMatMult(Btilde_pred_list[[modelidx]], t(thetatilde)) # out of sample point estimation.

      if(saveparams == TRUE){
        sigma2history[iter - nburnin] = sigma_2 
        lambdahistory[iter - nburnin] = lambda
        tauhistory[iter - nburnin] = tau
        Jhistory[iter - nburnin] = modelidx + 3 # first (the basic) model: sp{1,x,x^2,x^3}: J=4.
      }
    }
  } # end for loop (MCMC iteration)
  
  
  
  
  ### 5) posterior mean and credible interval of the estimaed function
  model_averaging = rowMeans(testpredictions) # posterior means of estimated out of sample function values.
  
  for(j in 1:n_pred){
    jth_testpredictions = testpredictions[j,]
    
    upper[j] = quantile(jth_testpredictions, probs = 0.975) # upper quantile of pointwise interval estimation of jth test point
    lower[j] = quantile(jth_testpredictions, probs = 0.025) # lower quantile of pointwise interval estimation of jth test point
  }

  
  if(plot_fit == TRUE){
    plot(x, y, cex = 0.5, font.main = 2, xlab ="", ylab = "", main = "")
    title("Data (dotted), posterior mean (blue), 95% coverage (grey)", adj = 0.02, line = -1)
    polygon(c(rev(x_pred), x_pred), c(rev(lower), upper), col = adjustcolor("grey", alpha.f=0.5) , border = NA)
    lines(x_pred, model_averaging, col = "blue", lwd = 1.5)
  }
  
  
  if(saveparams == TRUE){
    return(list("x" = x,
                "y" = y,
                "x_pred" = x_pred,
                "y_pred" = model_averaging,
                "upper" = upper,
                "lower" = lower,
                "J" = Jhistory,
                "sigma2" = sigma2history,
                "lambda" = lambdahistory,
                "tau" = tauhistory))
  }else{
    return(list("x" = x,
                "y" = y,
                "x_pred" = x_pred,
                "y_pred" = model_averaging,
                "upper" = upper,
                "lower" = lower))
  }
  
  
}
  