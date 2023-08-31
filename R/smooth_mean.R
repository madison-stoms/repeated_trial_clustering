
# input estimated mean and smooth using Gaussian kernel
smooth_mean = function(mu, K, L, t, nknots) {
  
  # loop over mu elements 
  musmooth = list()
  for(k in 1:K) {
    mutemp = matrix(NA, nrow =length(t), ncol = L)
    for(l in 1:L){
      
      temp = mu[[k]][,l]
      mutemp[,l] = smooth.spline(t, temp, nknots = nknots,  cv = TRUE)$y
      
    }
    musmooth[[k]] = mutemp
  }
  
  return(musmooth)
  
}
