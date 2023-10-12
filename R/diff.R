diff <-function(m,K){
  
# parameter  m: difference order
# parameter  K: size  

M = matrix(0,nrow=K-m,ncol=K)
c = rep(0,m+1)

for(i in 0:m)
c[i+1] = (-1)^(i+1)*factorial(m)/(factorial(i)*factorial(m-i))

for(i in 1:(K-m)) M[i,i:(i+m)] = c

return(M)
}
