# Stage 1 Estep: estimate tau
Estep = function(dat, theta) {
  
  # extract theta elements 
  mu = theta$mu; pi = theta$pi; 
  sigma2 = theta$sigma2; W = theta$W; lambda = theta$lambda
  K = length(pi)
  
  # extract elements
  m = length(unique(dat$trial)); n = length(unique(dat$subj)); 
  L = length(unique(dat$trialclus)); p = length(unique(dat$t))
  trialclus = sort(unique(trialclus))
  
  Ywide = dat %>% 
    pivot_wider(id_cols = c(subj, trial), names_from = t, values_from = y, names_prefix = "t")
  
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
        l = trialclus[j]
        
        # find c comps
        Yij = dat %>% filter(subj == i, trial == j) %>% pull(y)
        cvec[j] = t(Yij - mu[[k]][,l]) %*% Siginv[[k]][[l]] %*% (Yij - mu[[k]][,l])
        
      }
      C[i,k] = sum(cvec)
    }
  }
  
  
  
  # start responsibilty loop
  tau = matrix(NA, nrow = n, ncol = K)
  for(i in 1:n) {

    m = c()
    for(l in 1:L) { m[l] = sum(trialclus == l) }
    
    A = matrix(NA, nrow = K, ncol = K)
    for(q in 1:K) {
      for(s in 1:K) {
        
        bprod = 1
        for(l in 1:L) { brpod = B[[l]][q,s]^(-m[l]/2) * bprod }
        
        A[q,s] = bprod * exp(C[i,s]/2 - C[i,q]/2) * (pi[q]/pi[s])
        
      } }
    
    tau[i,] = 1 / colSums(A)
    
  }
  
  # check for missingness and reinitialize
  tau = ifelse(is.na(tau), runif(1), tau)
  tau = tau / rowSums(tau)
  
  return(list(tau = tau))
  
}

