# Stage 1 Estep: estimate tau
Estep = function(Ylong, t, theta) {
  
  # extract theta elements 
  mu = theta$mu; pi = theta$pi; 
  sigma2 = theta$sigma2; W = theta$W; lambda = theta$lambda
  K = length(pi)
  
  # extract elements
  m = length(unique(Ylong$trial)); n = length(unique(Ylong$subj)); 
  L = length(unique(Ylong$trialclus)); p = length(t)
  trialcluslong = Ylong$trialclus
  trialclus = trialcluslong[1:m]
  
  # extract just y vals
  Y = Ylong %>%
    dplyr::select(-c(subj, trialclus, trial)) %>% as.matrix()
  
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
  
  
  # find eta matrix for responsibility loop
  eta = matrix(NA, nrow = n, ncol = K)
  for(k in 1:K) {
    for(i in 1:n) {
      etavec = c()
      for(j in 1:m) {
        
        # identify trial type
        l = trialclus[j]
        
        # find c comps
        Yij = Ylong %>% filter(subj == i, trial == j) %>% dplyr::select(-c(subj,trialclus,trial)) %>% .[1,] %>% as.numeric()
        etavec[j] = matrix(Yij - mu[[k]][,l], nrow = 1) %*% Siginv[[k]][[l]] %*% matrix(Yij - mu[[k]][,l], ncol = 1)
        
      }
      eta[i,k] = sum(etavec)
    }
  }
  
  
  
  # start responsibility loop
  tau = matrix(NA, nrow = n, ncol = K)
  for(i in 1:n) { # start subj loop
    
    # extract how many trials of each type
    Lsum = c(); for(l in 1:L) { Lsum[l] = sum(trialclus == k) }
    
    for(k in 1:K) { # start cluster loop
      tauinner = 0
      for(c in 1:K) { # start inner cluster loop
        if(c == k) next # skip when sum = 0
        
        # calculate prod over trials 
        detprod = 1 
        for(l in 1:L) { detprod = detprod * (Sigdet[[c]][[l]] / Sigdet[[k]][[l]])^(-Lsum[l]/2) }
        
        # calculate 1/tau[i,k]
        tauinner = tauinner + ( pi[c]/pi[k] * exp((eta[i,k] - eta[i,c])/2) * detprod )
        
      }
      if(tauinner == 0) { tauinner = 0.0000001 }
      if(is.infinite(tauinner)) { tauinner = 9999999 }
      tau[i,k] = 1 / tauinner
    }
  }

  # check for missingness and reinitialize
  tau = ifelse(is.na(tau), runif(1), tau)
  tau = tau / rowSums(tau)
  
  
  return(list(tau = tau))
  
}





































