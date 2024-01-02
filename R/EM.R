
# EM Algorithm 
EM = function(dat, K, stoppt, cov.knots, mu.knots) {
  
  # arrange dat
  dat = dat %>%
    arrange(subj, trial, t)
  
  # extract elements
  n = length(unique(dat$subj))
  L = length(unique(dat$trial_type))
  t = unique(dat$t)
  
  # wide version of data
  Ywide = dat %>% 
    mutate(t = as.character(t)) %>%
    pivot_wider(id_cols = c(subj, trial_type, trial), names_from = t, values_from = y, names_prefix = "t") 
  
  # initialize tau
  #clus.initial = kmeans(Yclus, K)$cluster
  clus.initial = sample(1:K, n, replace = TRUE)
 # clus.initial = initial
  #clus.initial = rep(1:K, each = n/2)
  tau.initial = as.data.frame(model.matrix(~ -1 + factor(clus.initial)))
  # tau.initial = ifelse(tau.initial == 1, 0.7, 0.3)
  delta.initial = list(tau = tau.initial)
  
  # calculate rho
  typelong = Ywide$trial_type
  v = as.data.frame(model.matrix(~ -1 + factor(typelong)))
  rho = colSums(v) / nrow(v)
  
  # intial iteration
  theta = Mstep(Ywide, t, cov.knots, mu.knots, delta.initial, rho, r = 1)
  delta = Estep(Ywide, t, theta)
  deltaold = delta
  
  # first E and M steps
  for(i in 1:3) {

    theta = Mstep(Ywide, t, cov.knots, mu.knots, delta, rho, r = 1)
    delta = Estep(Ywide, t, theta)
    deltaold = delta
    #print(i + 1)

  }
  
  # start EM loop
  taulist = list()
  tol = 100; iter = 1
  while(tol > stoppt) {

    # E and M steps
    theta = Mstep(Ywide, t, cov.knots, mu.knots, delta, rho)
    delta = Estep(Ywide, t, theta)
    
     # taulist[[iter]] = list(tau = delta$tau, initial = clus.initial,
    #                        clus = max.col(delta$tau), theta = theta)
    # 
    # check convergence
    tol = max(abs(deltaold$tau - delta$tau)) 
    
    # prepare for next loop
    iter = iter + 1
    deltaold = delta
    
    if(iter == 20) { break }
    
  }
 # theta = Mstep(Ywide, t, cov.knots, mu.knots, delta, rho)
  
# BIC_check = funBIC(Ywide, theta, K = K)
  
  # add clus to Ywide
  Ywide = tibble(subj = unique(dat$subj), clus = max.col(delta$tau)) %>%
    full_join(Ywide, by = "subj")
  
  return(list(clus = max.col(delta$tau), delta = delta, theta = theta, # BIC = BIC_check,
              clus.initial = clus.initial, t = t, trial_type = typelong, iter = iter + 2,
              n = n, K = K, L = L, Ywide = Ywide))
  
  return(taulist)
  
}
