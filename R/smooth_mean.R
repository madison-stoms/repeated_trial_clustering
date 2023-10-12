
# input estimated mean and smooth using Gaussian kernel
smooth_mean = function(mukl, t, nknots) {
  
  musmooth = smooth.spline(t, mukl, nknots = nknots, cv = TRUE)$y
      
  return(musmooth)
  
}
