
# EM Algorithm 
EM = function(dat, K, stop, cov.knots, mu.knots) {
  
  # arrange dat
  dat = dat %>%
    arrange(subj, trial, t)
  
  # extract elements
  n = length(unique(dat$subj))
  L = length(unique(dat$trial_type))
  t = unique(dat$t)
  
  # wide version of data
  Ywide = dat %>% 
    pivot_wider(id_cols = c(subj, trial_type, trial), names_from = t, values_from = y, names_prefix = "t") 
  
  # initialize tau
  #clus.initial = kmeans(Yclus, K)$cluster
  clus.initial = sample(1:K, n, replace = TRUE)
  #clus.initial = rep(1:K, each = n/2)
  tau.initial = as.data.frame(model.matrix(~ -1 + factor(clus.initial)))
  # tau.initial = ifelse(tau.initial == 1, 0.7, 0.3)
  delta.initial = list(tau = tau.initial)
  
  # calculate rho
  typelong = Ywide$trial_type
  v = as.data.frame(model.matrix(~ -1 + factor(typelong)))
  rho = colSums(v) / nrow(v)
  
  # first E and M steps
  theta = Mstep(Ywide, t, cov.knots, mu.knots, delta.initial, rho)
  delta = Estep(Ywide, t, theta)
  deltaold = delta
  
  # start EM loop
  tol = 100; iter = 0
  while (tol > stop) {
    
    if(iter > 10) { break }

    # E and M steps
    theta = Mstep(Ywide, t, cov.knots, mu.knots, delta, rho)
    delta = Estep(Ywide, t, theta)
    
    # check convergence
    tol = max(abs(deltaold$tau - delta$tau)) 
    
    # prepare for next loop
    iter = iter + 1
    deltaold = delta
    
  }
  
 # BIC_check = funBIC(Y, theta, rho, K = K, L = L, n = n, m = m)
  
  # add clus to Ywide
  Ywide = tibble(subj = unique(dat$subj), clus = max.col(delta$tau)) %>%
    full_join(Ywide, by = "subj")
  
  return(list(clus = max.col(delta$tau), delta = delta, theta = theta, #BIC = BIC_check, 
              clus.initial = clus.initial, t = t, trial_type = typelong, iter = iter + 2,
              n = n, K = K, L = L, Ywide = Ywide))
  
}
