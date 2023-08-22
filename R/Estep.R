# Stage 1 Estep: estimate tau
Estep = function(Y, ID, trial, trialind, theta) {
  
  # extract elements 
  mu = theta$mu; pi = theta$pi; 
  sigma2 = theta$sigma2; W = theta$W; lambda = theta$lambda
  m = length(trialind); p = ncol(Y); n = length(unique(ID))
  K = length(pi); L = length(unique(trialind))
  
  
  # loop for variances and inverses
  Sig = list(); Siginv = list(); Sigdet = list()
  for(k in 1:K) {
    Sigtemp = list(); Siginvtemp = list(); Sigdettemp = list()
    for(l in 1:L) {
      
      r = ncol(W[[k]][[l]])
      
      if(length(lambda[[k]][[l]]) == 1) {
        diaglam = as.matrix(lambda[[k]][[l]])
      }
      else {
        diaglam = diag(lambda[[k]][[l]])
      }
      
      Sigtemp[[l]] = W[[k]][[l]] %*% diaglam %*% t(W[[k]][[l]]) + sigma2 * diag(p)
      Siginvtemp[[l]] = solve(Sigtemp[[l]])
      Sigdettemp[[l]] = det(sqrt(diaglam) %*% t(W[[k]][[l]]) %*% W[[k]][[l]] %*% sqrt(diaglam) + sigma2 * diag(r))
      
    }
    Sig[[k]] = Sigtemp; Siginv[[k]] = Siginvtemp; Sigdet[[k]] = Sigdettemp
    
  }
  
  
  
  # find B matrix
  B = list()
  for(l in 1:L) {
    
    bmat = matrix(NA, nrow = K, ncol = K)
    for(q in 1:K) { # dummy index
      for(s in 1:K) { # dummy index
        
        bmat[q,s] = Sigdet[[q]][[l]] / Sigdet[[s]][[l]]
        
      } }
    B[[l]] = bmat
    
  }
  
  
  # find C matrix
  C = matrix(NA, nrow = n, ncol = K)
  for(k in 1:K) {
    for(i in 1:n) {
      cvec = c()
      for(j in 1:m) {
        
        # identify trial type
        l = trialind[j]
        
        # find c comps
        Yij = Y[which(ID == i & trial == j),]
        cvec[j] = t(Yij - mu[[k]][,l]) %*% Siginv[[k]][[l]] %*% (Yij - mu[[k]][,l])
        
      }
      C[i,k] = sum(cvec)
    }
  }
  
  
  
  # start responsibilty loop
  tau = matrix(NA, nrow = n, ncol = K)
  for(i in 1:n) {
    m1 = sum(trialind == 1); m2 = sum(trialind == 2)
    
    A = matrix(NA, nrow = K, ncol = K)
    for(q in 1:K) {
      for(s in 1:K) {
        
        A[q,s] = B[[1]][q,s]^(-m1/2) * B[[2]][q,s]^(-m2/2) * exp(C[i,s]/2 - C[i,q]/2) * (pi[q]/pi[s])
        
      } }
    
    tau[i,] = 1 / colSums(A)
    
  }
  
  # check for missingness and reinitialize
  tau = ifelse(is.na(tau), runif(1), tau)
  tau = tau / rowSums(tau)
  
  return(list(tau = tau))
  
}

