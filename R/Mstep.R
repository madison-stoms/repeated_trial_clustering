
# Stage 1 Mstep: estimate pi and mu
Mstep = function(Y, ID, trial, trialind, t, cov.knots, delta) {
  
  # extract elements
  tau = delta$tau; m = length(trialind); p = ncol(Y); n = length(unique(ID))
  K = ncol(tau); L = length(unique(trialind))
  
  # calculate rho
  v = as.data.frame(model.matrix(~ -1 + factor(trialind)))
  vsum = colSums(v)
  rho = vsum / m
  
  
  # update pi
  tausum = colSums(tau)
  pi = tausum / n
  
  
  # update mu
  muraw = list()
  for(k in 1:K) {
    
    # loop over trials
    mutemp = matrix(NA, nrow = p, ncol = L)
    for(l in 1:L) {
      
      taulong = rep(tau[,k], each = m); vlong = rep(v[,l], n)
      mutemp[,l] = apply(Y * taulong * vlong, 2, sum) / (tausum[k] * vsum[l])
      
    }
    muraw[[k]] = mutemp
    
  }
  mu = smooth_mean(muraw, K = K, L = L, t = t)
  
  
  
  
  # estimate W and sigma2
  W = list(); lambda = list(); sigma2 = 0
  for(k in 1:K) {
    
    Wtemp = list(); lambdatemp = list()
    for(l in 1:L) {
      
      # calculate weighted covariance
      taulong = rep(tau[,k], each = m); vlong = rep(v[,l], n)
      Yc = sweep(Y, 2, mu[[k]][,l], "-")
      Sraw = (t(Yc * taulong * vlong) %*% Yc) / (pi[k] * n * rho[l] * m)
      
      # smooth covariance
      mod = fbps.cov(Sraw, knots = cov.knots)
      S = mod$cov
      sigma2kl = mod$var
      
      # decompose smoothed S
      decomp = eigen(S)
      evals = decomp$values
      efun = decomp$vectors
      
      var_prop =  cumsum(evals) / sum(evals)
      r = min(which(var_prop > 0.96))
      
      lambdatemp[[l]] = evals[1:r]
      Wtemp[[l]] = efun[,1:r] 
      
      # fix matrix problem
      if(is.null(dim(Wtemp[[l]]))) {
        Wtemp[[l]] = as.matrix(Wtemp[[l]], ncol = 1)}
      
      
      # weighted sum 
      sigma2 = sigma2kl * pi[k] * rho[l] + sigma2
      
    }
    W[[k]] = Wtemp; lambda[[k]] = lambdatemp
  }
  
  
  #return(efun)
  return(list(pi = pi, mu = mu, W = W, lambda = lambda, sigma2 = sigma2))
  
}


