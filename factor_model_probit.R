# fit a probit factor model  
scotch <- read.csv("scotch.csv")
subdat=subset(scotch,select = - Other.Brands)
oneto20=seq(1,20)
tofirst=c(11,9,16,18)
reor=c(tofirst,oneto20[-tofirst])
subdat=subdat[,reor]
Y=t(as.matrix(subdat))
Z=Y
Y=as.numeric(Y)


## function to update beta
update_beta <- function(X,Y,sigma,v = 1/100,beta0 = 0){
V0inv = diag(v,ncol(X),ncol(X))
a = 0.1; b = 0.1
p = ncol(X)
beta0 = rep(beta0,ncol(X))
XX = t(X)%*%X
beta_hat = solve(XX)%*%t(X)%*%Y
V1 = solve((1/sigma^2)*XX + V0inv)
beta1 = V1%*%((1/sigma^2)*XX%*%beta_hat + V0inv%*%beta0)

L = t(chol(V1))
beta = beta1 + L%*%rnorm(p)
return(beta)
}

n = ncol(Z)
p = nrow(Z)
k = 2 # number of factors

B = 0.1*matrix(rnorm(p*k),p,k)
diag(B) = abs(diag(B))
f = 0.1*matrix(rnorm(k*n),k,n)
mu = rep(0,p)
z = Z
mc = 1000

Bsave = array(0,c(p,k,mc))
mus = array(0,c(p,mc))
for (j in 1:(k-1)){
  B[j,(j+1):k] = 0
}
for (iter in 1:mc){
  for (j in 1:p){
    ind = 1:min(j,k)
    temp = update_beta(t(rbind(rep(1,n),f[ind,])),z[j,],1)
    B[j,ind] = temp[2:length(temp)]
    mu[j] = temp[1]
  }
  
  Ytemp = diag(1,p,p)%*%(z - mu)
  Xtemp = diag(1,p,p)%*%B
  
  for (i in 1:n){
    f[,i] = update_beta(Xtemp,Ytemp[,i],1,v = 1)
  }
  
  mutemp = B%*%f+mu
  a = rep(0,n*p)
  b = rep(1,n*p)
  a[Y==1] = pnorm(0,mutemp[Y==1],1)
  b[Y==0] = pnorm(0,mutemp[Y==0],1)
  u = runif(n*p,a,b)
  z = qnorm(u,mutemp,1)
  z = matrix(z,p,n)
  Bsave[,,iter] = B
  mus[,iter] = mu
}

Bmean = as.numeric(apply(Bsave[,,seq(200,mc)],c(1,2),mean))
Bmean = matrix(Bmean,p,k)
plot(Bmean[,1],Bmean[,2],pch='.')
cname=colnames(subdat)
text(Bmean[,1],Bmean[,2],cname,cex=.8)

