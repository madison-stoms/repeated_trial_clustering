
# calculate BIC
funBIC = function(Ywide, theta, K) {
  
  W = theta$W; mu = theta$mu; sigma2 = theta$sigma2
  p = nrow(W[[1]][[1]]); pi = theta$pi
  L = length(unique(Ywide$trial_type))
  
  # get sig stuff
  Sig = list(); r = c(); iter = 1
  for(k in 1:K) {
    
    Siginner = list()
    for(l in 1:L) {
      Siginner[[l]] =  W[[k]][[l]] %*% t(W[[k]][[l]]) + sigma2 * diag(p)
      r[iter] = ncol(W[[k]][[l]])
      iter = iter + 1
    }
    Sig[[k]] = Siginner
    
  }
  
  # loop for loglik
  loglik = 0
  for(i in 1:nrow(Ywide)) {
    
    Yijwide = Ywide[i,]; l = Yijwide$trial_type
    Yij = dplyr::select(Yijwide, -c(subj, clus, trial_type, trial)) %>%
      as.numeric()
    
    for(k in 1:K) {
      dclus = mvnfast::dmvn(Yij, mu[[k]][,l], Sig[[k]][[l]], log = TRUE)
      if(dclus == 0) { dclus = 0.000000001 }
      loglik = pi[k] * dclus + loglik
    }
  }
  
  # BIC
  ravg = round(mean(r))
  bic = -2 * loglik + log(nrow(Ywide)*p) * (2*K*ravg + 2*L*ravg -2*ravg + K + L - 1)
  
  return(bic)
      
}

 