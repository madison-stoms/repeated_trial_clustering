
reconstruct_curves = function(res) {
  
  # get scores
  scores = reconstruct_scores(res)
  
  # pull subjs 
  subjs = unique(res$Ywide$subj)
  
  # reconstruct curves
  recon = matrix(NA, nrow = nrow(res$Ywide), ncol = length(res$t)); iter = 1
  for(i in 1:res$n) {
    
    dati = filter(res$Ywide, subj == subjs[i])
    typei = dati$trial_type
    
    for(j in 1:nrow(dati)) {
      
      k = res$clus[i]; l = typei[j]
      recon[iter,] = res$theta$W[[k]][[l]] %*% scores[[k]][[l]][iter,] + res$theta$mu[[k]][,l]
      iter = iter + 1
      
    }
  }
  
  recon_dat = as_tibble(recon) %>%
    mutate(subj = res$Ywide$subj, trial = res$Ywide$trial,
           trial_type = res$Ywide$trial_type, clus = res$Ywide$clus) %>%
    pivot_longer(cols = starts_with("V"), names_to = "t", values_to = "y") %>%
    mutate(t = rep(res$t, nrow(res$Ywide))) %>%
    mutate(clus = factor(clus), trial_type = factor(trial_type)) %>%
    mutate(ID = str_c("Subj ", subj, ", Trial ", trial, ", Type ", trial_type))
  
  return(recon_dat)
}
