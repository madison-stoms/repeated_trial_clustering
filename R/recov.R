recov <-function(Sig,ratio){

## bin a covariance matrix while ignoring its diagonal elements

## Sig: an estimated (or true) covariance matrix 
## ratio: the bin size

m=dim(Sig)[1]
diag(Sig)=rep(NA,m)
m1=sqrt(length(Sig))/ratio
bin_cov=matrix(0,m1,m1)

if(ratio==1){
bin_cov = Sig
Diag = rep(0,m)
for(i in 2:(m-1))
  Diag[i] = (bin_cov[i-1,i]+bin_cov[i+1,i]+bin_cov[i,i-1]+bin_cov[i,i+1])/4
Diag[1] = (bin_cov[1,2]+bin_cov[2,1])/2
Diag[m] = (bin_cov[m-1,m]+bin_cov[m,m-1])/2
diag(bin_cov) = Diag
  }

if(ratio>=2){

for (i in 1:(m1-1))
for (j in (i+1):m1)
bin_cov[i,j]=mean(Sig[(i-1)*ratio+(1:ratio),(j-1)*ratio+(1:ratio)])

bin_cov=bin_cov+t(bin_cov)

for (i in 1:m1)
bin_cov[i,i]=mean(Sig[(i-1)*ratio+(1:ratio),(i-1)*ratio+(1:ratio)], na.rm=1)
}

return(bin_cov)
}
