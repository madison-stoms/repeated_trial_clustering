
# simulate clustered data using single PCs
simulate_data = function(n, m, t, pi = c(0.5,0.5), rho = c(0.5, 0.5), mu, W, lambda, sigma2) {
  
  # define indicies
  K = length(pi); L = length(rho); p = length(t)
  
  # construct true cluster assignments
  clus = c()
  for(k in 1:K ) { clus = c(clus, rep(k, n *pi[k])) }
  
  # construct trial indicators
  tind = c()
  for(l in 1:L ) { tind = c(tind, rep(l, m *rho[l])) } 
  
  
  # simulate data
  X = tibble(); iter = 1
  for(i in 1:n ) {
    
    for(j in 1:m ) {
      
      # define group
      ii = clus[i]
      jj = tind[j]
      
      # calculate curve
      temp = as.numeric(mu[[ii]][,jj] + W[[ii]][[jj]] * rnorm(1, 0, sqrt(lambda[[ii]][,jj])) + rnorm(p, 0, sqrt(sigma2)))
      X = rbind(X, tibble(ID = iter, subj = i, trial = j, subjclus = ii, trial_type = jj, t = t, y = temp))
      
      iter = iter + 1
      
    }
  }
  
  return(list(X = X, clus = clus))
  
}
