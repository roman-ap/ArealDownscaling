########## Spatial downscaling I
setwd("~/Nextcloud/UofE/Project/Belgium")

### Data sets
load("./BE_LAU2.RData")
load("./BE_LAU2pop.RData")
load("./BE_grid.RData")

### Packages
library(sf)
library(dplyr)
library(Matrix)

### Inputs
## Dimensions
N<-nrow(BE_grid)
M<-nrow(BE_LAU2)
## Overlap
BEgridded_LAU2<-st_intersection(BE_grid,BE_LAU2) %>%
  mutate(Area_km2=as.numeric(st_area(GEOMETRY))/1000000) %>%
  mutate(GC_ID=as.character(GC_ID),LAU2_ID=as.character(LAU2_ID)) %>%
  st_set_geometry(NULL)
## Matrix A
A_match<-st_relate(BE_LAU2,BE_grid,pattern="T********",sparse=TRUE)
Totals<-BEgridded_LAU2 %>%
  group_by(GC_ID) %>%
  summarise(Total=sum(Area_km2))
Areas<-left_join(BEgridded_LAU2,Totals,by="GC_ID") %>%
  mutate(Proportion=Area_km2/Total)
A<-sparseMatrix(i=c( rep(1:M,sapply(A_match,length)) ),
                j=c( unlist(A_match) ),
                x=c( Areas$Proportion ),
                dims=c(M,N))
## Aggregate counts
y<-Matrix(data=BE_LAU2pop$Population,nrow=M,ncol=1)

### Priors
## Spatial effects
# Mean vector (358.6969=national population density)
PD_grid<-mutate(BEgridded_LAU2,AW_est=Area_km2*358.6969) %>%
  group_by(GC_ID) %>%
  summarise(Est=sum(AW_est))
mu_x<-Matrix(PD_grid$Est,nrow=N,ncol=1)
# Precision matrix
W_match<-st_relate(BE_grid,BE_grid,pattern="F***1****",sparse=TRUE)
W<-sparseMatrix(i=c( 1:N,rep(1:N,sapply(W_match,length)) ),
                j=c( 1:N,unlist(W_match) ),
                x=c( sapply(W_match,length),rep(-1,length(unlist(W_match))) ),
                dims=c(N,N))
## Hyperparameters
#sigma_ref<-exp(sum(log(diag(solve(W,sparse=TRUE))))/(2*N))
a_k<-3
b_k<-2e+06

### Useful figures
Aty<-crossprod(A,y)
AtA<-crossprod(A)
v_x<-crossprod(W,mu_x)
q_x<-as.numeric(crossprod(mu_x,v_x))

approx_C<-function(x,logfx){
  n_c=length(x)
  h_c=x[2]-x[1]
  w_c=c(1,rep(c(4,2),(n_c-3)/2),4,1)
  return((h_c/3)*sum(w_c*exp(logfx)))
}

### Approximation of the posterior marginal of k
## Logarithm of the posterior marginal of k
logmarginalk<-function(k,t=1e-2){
  Q=forceSymmetric(k*W+t*AtA)
  L=Cholesky(Q)
  v=k*v_x+t*Aty
  u=solve(L,v,system="A")
  logn=dgamma(k,shape=a_k+(0.5*(N-1)),rate=b_k+0.5*q_x,log=T)
  logd=as.numeric(determinant(L,logarithm=T)$modulus)-0.5*as.numeric(crossprod(v,u))
  return(logn-logd-3704015000)
}
## Optimization 
optimize_k<-optim(par=5e-07,fn=logmarginalk,method="L-BFGS-B",
                  lower=1e-10,upper=1e-1,control=list(fnscale=-1,ndeps=1e-09),hessian=T)
k_mode<-optimize_k$par
log_k_mode<-logmarginalk(k_mode)
## Exploration 
zs<-Matrix(seq(-2.8,3.2,by=0.2))
# Final values
ks<-zs/sqrt(as.numeric(-optimize_k$hessian))+k_mode
logmarginalks<-sapply(ks,logmarginalk)
marginalk<-exp(logmarginalks-mean(logmarginalks))/approx_C(ks,logmarginalks-mean(logmarginalks))
k_data<-tibble(k=as.vector(ks),d_k=marginalk)
library(ggplot2)
library(latex2exp)
ggplot(k_data,aes(x=as.vector(ks),y=marginalk))+
  geom_point(color="red")+
  expand_limits(x=c(1e-07,2.5e-07))+
  labs(x=TeX("$\\kappa$"),y=TeX("$p(\\kappa | \\mathbf{y})$"))+
  theme_bw()
# Weights
k_weights1<-(ks[2]-ks[1])*marginalk
save(k_weights1,file="k_weights1.RData")


### Posterior marginal of x
## Means
library(doParallel)
cl<-makeCluster(8)
registerDoParallel(cl)
x_means1<-foreach(k=as.vector(ks),.combine='cbind',.packages=c('Matrix')) %dopar% {
  t=1e-02
  Q=forceSymmetric(k*W+t*AtA)
  L=Cholesky(Q)
  v=k*v_x+t*Aty
  u=solve(L,v,system="A")
  return( u )
}
stopCluster(cl)
# RData
save(x_means1,file="x_means1.RData")

## Samples
library(doRNG)
S<-5000
cl<-makeCluster(4)
registerDoRNG(27)
registerDoParallel(cl)
x_vars1<-foreach(k=as.vector(ks),.combine='cbind',.packages=c('Matrix')) %dopar% {
  t=1e-02
  Q=forceSymmetric(k*W+t*AtA)
  L=Cholesky(Q,LDL=F)
  z=Matrix(rnorm(n=N*S),nrow=N,ncol=S)
  s=solve(L,z,system="Lt")
  vars=apply(s,1,var)
  return( solve(L,vars,system="Pt") )
}
stopCluster(cl)
# RData
save(x_vars1,file="x_vars1.RData")
