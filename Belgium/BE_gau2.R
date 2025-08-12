########## Spatial downscaling II
setwd("~/Nextcloud/UofE/Project/Belgium")

### Data sets
load("./BE_LAU2.RData")
load("./BE_LAU2pop.RData")
load("./BE_grid.RData")
load("./BE.RData")

### Packages
library(sf)
library(dplyr)
library(Matrix)

### Inputs
## Dimensions
N<-nrow(BE_grid)
M<-nrow(BE_LAU2)
## Ancillary data
# CLC areas at grid cell level
myBEgridded<-st_set_geometry(BE,NULL) %>% 
  mutate(myCLC_class=recode(CLC_class,
                            '111' = 'Urban',
                            '112' = 'Rural',
                            '121' = 'Commercial',
                            '122' = 'Transport',
                            '123' = 'Transport',
                            '124' = 'Transport',
                            '131' = 'Sites',
                            '132' = 'Sites',
                            '133' = 'Sites',
                            '141' = 'Leisure',
                            '142' = 'Leisure',
                            '211' = 'Arable',
                            '212' = 'Arable',
                            '213' = 'Arable',
                            '221' = 'Groves',
                            '222' = 'Groves',
                            '223' = 'Groves',
                            '231' = 'Pasture',
                            '241' = 'Crops',
                            '242' = 'Crops',
                            '243' = 'Crops',
                            '244' = 'Crops',
                            '311' = 'Forest',
                            '312' = 'Forest',
                            '313' = 'Forest',
                            '321' = 'Vegetation',
                            '322' = 'Vegetation',
                            '323' = 'Vegetation',
                            '324' = 'Vegetation',
                            '331' = 'Bare',
                            '332' = 'Bare',
                            '333' = 'Bare',
                            '334' = 'Bare',
                            '335' = 'Bare',
                            '411' = 'Wetlands',
                            '412' = 'Wetlands',
                            '421' = 'Wetlands',
                            '422' = 'Wetlands',
                            '423' = 'Wetlands',
                            '511' = 'Water',
                            '512' = 'Water',
                            '521' = 'Water',
                            '522' = 'Water',
                            '523' = 'Water',
                            .default='Other')) %>%
  group_by(GC_ID,myCLC_class) %>%
  summarise(Area=sum(Area_km2)) %>%
  tidyr::spread(key=myCLC_class,value=Area,fill=0) %>%
  mutate(Inhabited=Urban+Rural,
         Rest=Arable+Groves+Pasture+Crops+Forest) %>%
  select(GC_ID,Inhabited,Urban,Rural,Rest)
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
# Mean vector (2110.596=national population density)
mu_x<-Matrix(2110.596*myBEgridded$Inhabited,nrow=N,ncol=1,sparse=TRUE)
# Uninhabited gridcells
zeros<-if_else(myBEgridded$Inhabited>0,F,T)
# Precision matrix
grid_match<-st_relate(BE_grid,BE_grid,pattern="F***1****",sparse=TRUE)
W_match<-grid_match[1:N]
W_match[which(zeros==TRUE)]<-lapply(grid_match[zeros==TRUE],setdiff,y=which(zeros==TRUE)) 
W<-sparseMatrix(i=c( which(zeros==TRUE), which(zeros==FALSE), rep(1:N,sapply(W_match,length)) ),
                j=c( which(zeros==TRUE), which(zeros==FALSE), unlist(W_match) ),
                x=c( rep(1e+10,sum(zeros)), sapply(W_match[which(zeros==FALSE)],length), rep(-1,length(unlist(W_match))) ),
                dims=c(N,N))
## Hyperparameters
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
  return(logn-logd-3703967850)
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
k_weights2<-(ks[2]-ks[1])*marginalk
save(k_weights2,file="k_weights2.RData")


### Posterior marginal of x
## Means
library(doParallel)
cl<-makeCluster(8)
registerDoParallel(cl)
x_means2<-foreach(k=as.vector(ks),.combine='cbind',.packages=c('Matrix')) %dopar% {
  t=1e-02
  Q=forceSymmetric(k*W+t*AtA)
  L=Cholesky(Q)
  v=k*v_x+t*Aty
  u=solve(L,v,system="A")
  return( u )
}
stopCluster(cl)
# RData
save(x_means2,file="x_means2.RData")

## Samples
library(doRNG)
S<-5000
cl<-makeCluster(4)
registerDoRNG(27)
registerDoParallel(cl)
x_vars2<-foreach(k=as.vector(ks),.combine='cbind',.packages=c('Matrix')) %dopar% {
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
save(x_vars2,file="x_vars2.RData")
