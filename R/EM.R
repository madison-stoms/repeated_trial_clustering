
# EM Algorithm 
EM = function(Y, ID, trial, trialind, K, stop, t, cov.knots) {
  
  # extract elements
  m = length(unique(trial)); n = length(unique(ID)); L = length(unique(trialind)); p = length(t)
  trialind = trialind[1:m]
  
  # initialize delta with first trial of every observation
  
  Yclus = as_tibble(Y) %>%
    mutate(ID = ID, trial = trial) %>%
    pivot_longer(cols = -c(ID, trial), names_to = "time", values_to = "val") %>%
    pivot_wider(id_cols = c(ID, trial, time), names_from = ID, values_from = val, names_prefix = "ID") %>%
    dplyr::select(-c(trial, time)) %>%
    as.matrix() %>% t()
  
  clus.initial = kmeans(Yclus, K)$cluster
  tau.initial = as.data.frame(model.matrix(~ -1 + factor(clus.initial)))
  delta.initial = list(tau = tau.initial)
  
  # define initial clusters and format
  clus = max.col(delta.initial$tau)
  cluslong = rep(clus, each = m)
  trialindlong = rep(trialind, n)
  
  # calculate rho
  v = as.data.frame(model.matrix(~ -1 + factor(trialind)))
  vsum = colSums(v)
  rho = vsum / m
  
  # first two E and M steps
  theta = Mstep(Y, ID, trial, trialind, t, cov.knots, delta.initial)
  delta = Estep(Y, ID, trial, trialind, theta)
  
  theta = Mstep(Y, ID, trial, trialind, t, cov.knots, delta)
  delta = Estep(Y, ID, trial, trialind, theta)
  deltaold = delta
  
  # start EM loop
  tol = 100; iter = 0
  while (tol > stop & iter < 10) {
    
    # E and M steps
    theta = Mstep(Y, ID, trial, trialind, t, cov.knots, delta)
    delta = Estep(Y, ID, trial, trialind, theta)
    
    # check convergence
    tol = max(abs(deltaold$tau - delta$tau)) 
    
    # prepare for next loop
    iter = iter + 1
    deltaold = delta
    
  }
  
  BIC_check = funBIC(Y, theta, rho, K = K, L = L, n = n, m = m)
  
  return(list(delta = delta, theta = theta, BIC = BIC_check, clus.initial = clus.initial, iter = iter + 2))
  
}
