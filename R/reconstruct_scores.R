

reconstruct_scores = function(res) {
  
  Y = res$Ywide %>% dplyr::select(-c(subj, trial, trial_type, clus))
  
  # reconstruct full scores
  scores.clus = list()
  for(k in 1:res$K) {
    
    scoresk = list()
    for(l in 1:res$L) {
      
      Yc = t(sweep(Y, 2, res$theta$mu[[k]][,l], "-"))
      scoresk[[l]] = t(t(res$theta$W[[k]][[l]]) %*% Yc)
      
    }
    scores.clus[[k]] = scoresk
    
  }

  
  # # extract based on group assignment 
  # scoreslist = list(); iter = 1
  # for(i in 1:res$n) {
  #   
  #   Ji = filter(res$Ywide, subj == 18, trial_type == 1) %>% nrow()
  #   scoreslist[[i]] = matrix(NA, nrow = res$n, ncol = Ji)
  #   
  #   for(j in 1:Ji) {
  # 
  #     scores[iter,] = scores.clus[[res$clus[i]]][[res$trialclus[j]]][iter,]
  #     iter = iter + 1
  # 
  #   }
  # }
  
  return(scores.clus)

}
