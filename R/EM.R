
# EM Algorithm 
EM = function(dat, K, stop, cov.knots, mu.knots) {
  
  dat = dat %>%
    arrange(subj, trial, t)
  
  # extract elements
  m = length(unique(dat$trial)); n = length(unique(dat$subj)); 
  L = length(unique(dat$trialclus)); t = unique(dat$t); p = length(t)
  
  # format Y
  Ylong = dat %>% 
    pivot_wider(id_cols = c(subj, trialclus, trial), names_from = t, values_from = y, names_prefix = "t") 
  
  Y = Ylong %>% dplyr::select(-c(subj, trialclus, trial)) %>% as.matrix()
  
  trialclus = Ylong %>% pull(trialclus) %>% .[1:m]
  trialcluslong = Ylong$trialclus
  
  #clus.initial = kmeans(Yclus, K)$cluster
  clus.initial = sample(1:K, n, replace = TRUE)
  #clus.initial = rep(1:K, each = n/2)
  tau.initial = as.data.frame(model.matrix(~ -1 + factor(clus.initial)))
 # tau.initial = ifelse(tau.initial == 1, 0.7, 0.3)
  delta.initial = list(tau = tau.initial)
  
  # calculate rho
  v = as.data.frame(model.matrix(~ -1 + factor(trialclus)))
  vsum = colSums(v)
  rho = vsum / m
  
  # first two E and M steps
  theta = Mstep(Ylong, t, cov.knots, mu.knots, delta.initial)
 # if(!is.list(theta)) { return(theta) }
  delta = Estep(Ylong, t, theta)
  
  theta = Mstep(Ylong, t, cov.knots, mu.knots, delta)
  delta = Estep(Ylong, t, theta)
  deltaold = delta
  
  # start EM loop
  tol = 100; iter = 0
  while (tol > stop) {
    
    if(iter > 10) { break }

    # E and M steps
    theta = Mstep(Ylong, t, cov.knots, mu.knots, delta)
    delta = Estep(Ylong, t, theta)
    
    # check convergence
    tol = max(abs(deltaold$tau - delta$tau)) 
    
    # prepare for next loop
    iter = iter + 1
    deltaold = delta
    
   # print(iter)
    
  }
  
  BIC_check = funBIC(Y, theta, rho, K = K, L = L, n = n, m = m)
  
  return(list(clus = max.col(delta$tau), delta = delta, theta = theta, BIC = BIC_check, 
              clus.initial = clus.initial, t = t, trialclus = trialclus, iter = iter + 2,
              m = m, n = n, K = K, L = L, Y = Y, Ylong = Ylong))
  
}
