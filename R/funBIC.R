
# caculate BIC
funBIC = function(Y, theta, rho, K, L, n, m) {
  
  W = theta$W; mu = theta$mu; sigma2 = theta$sigma2; pi = theta$pi
  p = nrow(W[[1]][[1]])
  
  outerlik = c(); iter = 0; r = c()
  for(i in 1:n) {
    for(j in 1:m) {
      
      innerlik = 0; iter = iter + 1
      for(k in 1:K) {
        for(l in 1:L) {
          
          r[iter] = ncol(W[[k]][[l]])
          Yij = Y[iter,]
          
          sig =  W[[k]][[l]] %*% t(W[[k]][[l]]) + sigma2 * diag(p)
          dclus = mvnfast::dmvn(Yij, mu[[k]][,l], sig, log = FALSE)
          
          innerlik = pi[k] * rho[l] * dclus + innerlik 
          
        }
      }
      
      outerlik[i] = log(innerlik)
      
    }
  }
  loglik = sum(outerlik)
  
  
  # BIC
  ravg = round(mean(r))
  bic = -2 * loglik + log(nrow(Y)*p) * (2*K*ravg + 2*L*ravg -2*ravg + K + L - 1)
}
