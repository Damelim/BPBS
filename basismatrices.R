library(splines)

get_B_matrix = function(x, x_pred, degree = 3, num_interior_knots){ 
  
  # span(1, B_2 - C_2, ..., B_m - C_m)
  
  knots = seq(0, 1, length = num_interior_knots + 2)
  knots = knots[-1] ; knots = knots[-length(knots)]
  
  B = bs(x, knots = knots, degree = degree, intercept= T)
  
  colmeanvec = colMeans(B)
  n = nrow(B) ; p = ncol(B)
  colmeanmat = matrix(rep(colmeanvec, n), n, p, byrow = T)
  B = B - colmeanmat
  B[,1] = 1
  
  if(is.null(x_pred) == FALSE){
    n_pred = length(x_pred)
    B_pred = bs(x_pred, knots = knots, degree = degree, intercept= T)
    colmeanmat = matrix(rep(colmeanvec, n_pred), n_pred, p, byrow = T)
    B_pred = B_pred - colmeanmat
    B_pred[,1] = 1
  }else{
    B_pred = NULL
  }
  
  return(list("B" = B,  "B_pred" = B_pred))
  
}

get_P_matrix = function(n_col, deriv){  # P = D^T D and then, fill the first row and column as 0. 
  D = as.matrix(diff(diag(1, n_col), diff = deriv)) 
  P = crossprod(D) #eigenMapMatMult(t(D), D)
  P[1,1:ncol(P)] = 0 
  P[1:nrow(P),1] = 0
  return(P)
}