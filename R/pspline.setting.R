pspline.setting <-function(x,knots=35,p=3,m=2){
  
# x: the marginal data points
# K: the list of knots or the numbers of knots
# p: degrees for B-splines, with defaults values 3
# m: orders of difference penalty, with default values 2
#library(splines)
if(length(knots)==1)
{
  K = knots
  knots=seq(-p,K+p,length=K+1+2*p)/K
  knots = knots*(max(x)-min(x)) + min(x)
  }

if(length(knots)>1) 
  { knots = knots
    K = length(knots)-2*p-1
}


P = diff(m,K+p)
P = t(P)%*%P

### design matrix and some pre-calculation 
### for the penalty without interaction
### The idea about pre-calculation, when lambda
## is changed, only eigenvalues change.

B = spline.des(knots=knots, x=x, ord = p+1,outer.ok = TRUE)$design

Sig = t(B)%*%B
eSig = eigen(Sig)
V = eSig$vectors
E = eSig$values

Sigi_sqrt =V%*%diag(1/sqrt(E))%*%t(V)
Sigi = V%*%diag(1/E)%*%t(V)

tUPU = t(Sigi_sqrt)%*%P%*%Sigi_sqrt

Esig = eigen(tUPU)
U = Esig$vectors
s = Esig$values
s[(K+p-m+1):(K+p)]=0
A = B%*%Sigi_sqrt%*%U

List = list(
        "A" = A,
        "s" = s,
        "Sigi.sqrt"=Sigi_sqrt,
        "U" = U)

return(List)
}
