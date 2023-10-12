
# Stage 1 Mstep: estimate pi and mu
Mstep = function(Ywide, t, cov.knots, mu.knots, delta, rho) {
  
  # extract delta elements
  tau = delta$tau; K = ncol(tau)
  
  # extract elements
  n = length(unique(Ywide$subj)); p = length(t)
  L = length(unique(Ywide$trial_type))
  subjs = unique(Ywide$subj)
  
  # extract just y vals
  Y = Ywide %>%
    dplyr::select(-c(subj, trial_type, trial)) %>% as.matrix()
  
  # update pi
  tausum = colSums(tau)
  pi = tausum / n
  for(k in 1:K) { ifelse(pi[k] <= 3e-23, 0.01, pi[k]) }
  pi = pi / sum(pi)
  
  # calculate trial lengths
  Jvec = Ywide %>% 
    group_by(subj) %>%
    summarise(J = n()) %>% ungroup() %>% pull(J)
  
  # get trial type counts by subj
  typesums = Ywide %>% 
    group_by(subj, trial_type) %>%
    summarise(n = n()) %>% ungroup()
  
  # calculate v for vlongs
  v = as.data.frame(model.matrix(~ -1 + factor(Ywide$trial_type)))
  
  
  # calculate denominators 
  denom = list()
  for(k in 1:K) {
    denom_inner = list()
    for(l in 1:L) {
      inner = 0
      for(i in 1:44) { 
        typesumi = filter(typesums, subj == subjs[i], trial_type == l) %>% pull(n)
        inner = inner + tau[i,k] * typesumi
        if(is.nan(inner)) {stop(i)}
      }
      denom_inner[[l]] = inner
    }
    denom[[k]] = denom_inner
  }
  
  
  
  # # update mu
  # muraw = list()
  # for(k in 1:K) {
  #   
  #   taulong = c()
  #   for(i in 1:n) {taulong = c(taulong, rep(tau[i,k], Jvec[i]))}
  #   
  #   # loop over trials
  #   mutemp = matrix(NA, nrow = p, ncol = L)
  #   for(l in 1:L) {
  #     
  #     vlong = v[,l]
  #     mutemp[,l] = apply(Y * taulong * vlong, 2, sum) / denom[[k]][[l]]
  #     
  #   }
  #   muraw[[k]] = mutemp
  #   
  # }
  # mu = smooth_mean(muraw, K = K, L = L, t = t, nknots = mu.knots)
  # 
  # 
  
  
  # estimate mu, W, and sigma2
  mu = list(); W = list(); lambda = list(); sigma2 = 0
  for(k in 1:K) {
    
    taulong = c()
    for(i in 1:n) {taulong = c(taulong, rep(tau[i,k], Jvec[i]))}
    
    mutemp = matrix(NA, nrow = p, ncol = L); Wtemp = list(); lambdatemp = list()
    for(l in 1:L) {
      
      vlong = v[,l]
      
      # update mean
      temp = apply(Y * taulong * vlong, 2, sum) / denom[[k]][[l]]
      mutemp[,l] = smooth.spline(t, temp, nknots = mu.knots, cv = TRUE)$y
      
      # calculate weighted covariance
      Yc = sweep(Y, 2, mutemp[,l], "-")
      Sraw = (t(Yc * taulong * vlong) %*% Yc) / denom[[k]][[l]]
      
      # smooth covariance
      mod = fbps.cov(Sraw, knots = cov.knots)
      S = mod$cov
      S = (S + t(S)) / 2
      sigma2kl = mod$var
      
      # decompose smoothed S
      decomp = eigen(S)
      evals = decomp$values
      efun = decomp$vectors

      # choose trunc
      var_prop =  cumsum(evals) / sum(evals)
      r = min(which(var_prop > 0.97))
      
      # save results 
      lambdatemp[[l]] = evals[1:r]
      Wtemp[[l]] = efun[,1:r] 
      
      # fix matrix problem
      if(is.null(dim(Wtemp[[l]]))) {
        Wtemp[[l]] = as.matrix(Wtemp[[l]], ncol = 1)}
      
      # update weighted sum 
      sigma2 = sigma2kl * pi[k] * rho[l] + sigma2
      
    }
    
    mu[[k]] = mutemp; W[[k]] = Wtemp; lambda[[k]] = lambdatemp
    
  }
  
  return(list(pi = pi, mu = mu, W = W, lambda = lambda, sigma2 = sigma2))
  
}









