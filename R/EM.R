
# EM Algorithm 
EM = function(dat, K, stop, cov.knots) {
  
  # extract elements
  m = length(unique(dat$trial)); n = length(unique(dat$subj)); 
  L = length(unique(dat$trialclus)); p = length(unique(dat$t))
  trialclus = sort(unique(dat$trialclus))
  
  # initialize delta with first trial of every observation
  Yclus = dat %>%
    pivot_wider(id_cols = c(trial, t), names_from = subj, values_from = y, names_prefix = "subj") %>%
    dplyr::select(-c(trial, t)) %>%
    as.matrix() %>% t()
  
  clus.initial = kmeans(Yclus, K)$cluster
  tau.initial = as.data.frame(model.matrix(~ -1 + factor(clus.initial)))
  delta.initial = list(tau = tau.initial)
  
  # define initial clusters and format
  clus = max.col(delta.initial$tau)
  cluslong = rep(clus, each = m)
  trialcluslong = rep(trialclus, n)
  
  # calculate rho
  v = as.data.frame(model.matrix(~ -1 + factor(trialclus)))
  vsum = colSums(v)
  rho = vsum / m
  
  # first two E and M steps
  theta = Mstep(dat, cov.knots, delta.initial)
  delta = Estep(dat, theta)
  
  theta = Mstep(dat, cov.knots, delta)
  delta = Estep(dat, theta)
  deltaold = delta
  
  # start EM loop
  tol = 100; iter = 0
  while (tol > stop & iter < 10) {
    
    # E and M steps
    theta = Mstep(dat, cov.knots, delta)
    delta = Estep(dat, theta)
    
    # check convergence
    tol = max(abs(deltaold$tau - delta$tau)) 
    
    # prepare for next loop
    iter = iter + 1
    deltaold = delta
    
  }
  
  BIC_check = funBIC(Y, theta, rho, K = K, L = L, n = n, m = m)
  
  return(list(delta = delta, theta = theta, BIC = BIC_check, clus.initial = clus.initial, iter = iter + 2))
  
}
