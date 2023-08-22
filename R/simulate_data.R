
# simulate clustered data
simulate_data = function(n.t = n, m.t = m, K.t = K, L.t = L, r.t = r, t.t = t, pi = c(0.5,0.5), rho = c(0.5, 0.5), mu, W, lambda, sigma2) {
  
  # define observation interval length
  p = length(t.t)
  
  # construct true cluster assignments
  clus = c()
  for(k in 1:K.t) { clus = c(clus, rep(k, n.t*pi[k])) }
  
  # construct trial indicators
  tind = c()
  for(l in 1:L.t) { tind = c(tind, rep(l, m.t*rho[l])) } 
  
  
  # simulate data
  X = matrix(NA, nrow = n.t*m.t, ncol = (p + 3)); iter = 1
  for(i in 1:n.t) {
    
    eta = rnorm(r.t, 0, 9)
    for(j in 1:m.t) {
      
      # define group
      ii = clus[i]
      jj = tind[j]
      
      # calculate curve
      temp = mu[[ii]][,jj] + W[[ii]][[jj]] %*% matrix(mvrnorm(1, eta, diag(sqrt(lambda.t[[ii]][,jj]))), nrow = r.t) + rnorm(p, 0, sqrt(sigma2))
      X[iter,] = cbind(t(c(i, j, jj)), t(temp))
      
      iter = iter + 1
      
    }
  }
  
  return(list(X = X, clus = clus))
  
}