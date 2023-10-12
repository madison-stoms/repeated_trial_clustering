
# Stage 1 Estep: estimate tau
Estep = function(Ywide, t, theta) {
  
  # extract theta elements 
  mu = theta$mu; pi = theta$pi; 
  sigma2 = theta$sigma2; W = theta$W; lambda = theta$lambda
  K = length(pi); L = length(unique(Ywide$trial_type))
  
  # extract elements
  n = length(unique(Ywide$subj)); p = length(t)
  typelong = Ywide$trial_type
  subjs = unique(Ywide$subj)
  
  # extract just y vals
  Y = Ywide %>%
    dplyr::select(-c(subj, trial_type, trial)) %>% as.matrix()
  
  
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
      
      dati = filter(Ywide, subj == subjs[i]) 
      J = length(dati$trial)
      typei = dati$trial_type
      
      etavec = c()
      for(j in 1:J) {
        # identify trial type
        l = typei[j]
        # find c comps
        Yij = dati[j,] %>% dplyr::select(-c(subj, trial_type, trial)) %>% as.numeric()
        etavec[j] = matrix(Yij - mu[[k]][,l], nrow = 1) %*% Siginv[[k]][[l]] %*% matrix(Yij - mu[[k]][,l], ncol = 1)
      }
      eta[i,k] = sum(etavec)
    }
  }
 
  
  # start responsibility loop
  tau = matrix(NA, nrow = n, ncol = K)
  for(i in 1:n) { # start subj loop
    
    # extract how many trials of each type for subj i
    typei = filter(Ywide, subj == subjs[i]) %>% pull(trial_type)
    Lsum = c(); for(l in 1:L) { Lsum[l] = sum(typei == l) }
    
    for(k in 1:K) { # start cluster loop
      tauinner = 0
      for(c in 1:K) { # start inner cluster loop
        if(c == k) next # skip when sum = 0
        
        # calculate prod over trials 
        detprod = 1 
        for(l in 1:L) { detprod = detprod * ((Sigdet[[c]][[l]] / Sigdet[[k]][[l]])^(-Lsum[l]/2)) }
        if(detprod == 0) { detprod = 0.00001 }
        if(is.infinite(detprod)) { detprod = 999999999 }
        
        # calculate 1/tau[i,k]
        
        tauinner = tauinner + ( pi[c]/pi[k] * exp((eta[i,k] - eta[i,c])/2) * detprod )
        
        if(is.na(tauinner)) {print(c(i,k,c))}
      }
      if(tauinner < 0.0000001) { tauinner = 0.0000001 }
      if(tauinner > 999999999) { tauinner = 999999999 }
      tau[i,k] = 1 / tauinner
    }
  }

  # check for missingness and reinitialize
  tau = ifelse(is.na(tau), runif(1), tau)
  tau = tau / rowSums(tau)
  
  return(list(tau = tau))
  
}














