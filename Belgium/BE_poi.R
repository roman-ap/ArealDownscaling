########## Spatial downscaling via Laplace Approximations
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
library(matrixStats)


### Inputs
## Dimensions
N<-nrow(BE_grid)
M<-nrow(BE_LAU2)
P<-1L
## Ancillary data
# CLC areas at grid cell level
myBEgridded<-st_set_geometry(BE,NULL) %>% 
  mutate(myCLC_class=recode(CLC_class,
                            '111' = 'Urban',
                            '112' = 'Urban',
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
  mutate(Rest=Arable+Groves+Pasture+Crops+Forest) %>%
  select(GC_ID,Urban,Rest)
# Matrix Z
Z<-sparseMatrix(i=c( which(myBEgridded$Urban>0) ),
                j=c( rep(1,length(which(myBEgridded$Urban>0))) ),
                x=c( myBEgridded$Urban[which(myBEgridded$Urban>0)] ),
                dims=c(N,P))
## Overlap
BEgridded_LAU2<-st_intersection(BE_grid,BE_LAU2) %>%
  mutate(GC_ID=as.character(GC_ID),LAU2_ID=as.character(LAU2_ID)) %>%
  mutate(Area_km2=as.numeric(st_area(GEOMETRY))/1000000) %>%
  st_set_geometry(NULL)
## Matrix A
A_match<-st_relate(BE_LAU2,BE_grid,pattern="T********",sparse=TRUE)
A<-sparseMatrix(i=c( rep(1:M,sapply(A_match,length)) ),
                j=c( unlist(A_match) ),
                x=c( BEgridded_LAU2$Area_km2 ),
                dims=c(M,N))
## Aggregate counts
y<-Matrix(data=BE_LAU2pop$Population,nrow=M,ncol=1)

### Priors
## Fixed effects
mu_beta<-Matrix(10,P,1)
Q_beta<-Diagonal(P,1e-2)
## Spatial effects
# Mean vector
mu_u<-Matrix(0,nrow=N,ncol=1)
# Precision matrix
W_match<-st_relate(BE_grid,BE_grid,pattern="F***1****",sparse=TRUE)
W<-sparseMatrix(i=c( 1:N,rep(1:N,sapply(W_match,length)) ),
                j=c( 1:N,unlist(W_match) ),
                x=c( sapply(W_match,length),rep(-1,length(unlist(W_match))) ),
                dims=c(N,N))
## Hyperparameters
a_k<-2
b_k<-1e+2
## Useful figures
mu_x<-rbind2(mu_beta,mu_u)
q_u<-as.numeric(crossprod(mu_u,crossprod(W,mu_u)))

### Functions
eval_f<-function(x,k){
  Q_x=bdiag(Q_beta,k*W)
  v_x=crossprod(Q_x,x)
  beta=Matrix(x[1:P],1,P)
  u=Matrix(x[-c(1:P)],1,N)
  lambda=exp(tcrossprod(beta,Z)+u)
  rates=tcrossprod(A,lambda)
  f_xk=crossprod(y,log(rates))-sum(lambda)-0.5*crossprod(x,v_x)+crossprod(mu_x,v_x)
  return(as.numeric(f_xk))
}

eval_gradientf<-function(x,k){
  Q_x=bdiag(Q_beta,k*W)
  beta=Matrix(x[1:P],1,P)
  u=Matrix(x[-c(1:P)],1,N)
  lambda=exp(tcrossprod(beta,Z)+u)
  rates=tcrossprod(A,lambda)
  dlambda_dbeta=Z*as.vector(lambda)
  dh_dbeta=A%*%dlambda_dbeta
  dg_dbeta=crossprod(dh_dbeta,y/rates)-Matrix(colSums(dlambda_dbeta))
  dh_du=t(t(A)*as.vector(lambda))
  dg_du=crossprod(dh_du,y/rates)-Matrix(colSums(lambda))
  nabla_f=rbind2(dg_dbeta,dg_du)-crossprod(Q_x,x-mu_x)
  return(as.vector(nabla_f))
}

nonzeros<-function(column_number){
  regions=which(A[,column_number]!=0)
  if(length(regions)>1){
    block=colSums(A[regions,],sparseResult=TRUE)
    gridcells=which(block!=0)
  } else {gridcells=which(A[regions,]!=0)}
  return(union(column_number,gridcells))
}
h_nonzeros<-function(cells,dh_du,rates){
  d2h_nonzeros=dh_du[,cells]*as.vector(dh_du[,cells[1]])
  d2h_du2=crossprod(d2h_nonzeros,1/rates)
  return(as.vector(d2h_du2))
}
u_nonzeros<-apply(array(1:N),1,nonzeros)

eval_hessianf<-function(x,k){
  Q_x=bdiag(Q_beta,k*W)
  beta=Matrix(x[1:P],1,P)
  u=Matrix(x[-c(1:P)],1,N)
  lambda=exp(tcrossprod(beta,Z)+u)
  rates=tcrossprod(A,lambda)
  dlambda_dbeta=Z*as.vector(lambda)
  dh_dbeta=A%*%dlambda_dbeta
  dh_du=t(t(A)*as.vector(lambda))
  H_beta=Matrix(0,P,P)
  H_ubeta=Matrix(0,N,P)
  for (i in 1:P){
    d2lambda_dbeta2_i=dlambda_dbeta*as.vector(Z[,i])
    d2h_dbeta2_i=(A%*%d2lambda_dbeta2_i)*as.vector(rates)-dh_dbeta*as.vector(dh_dbeta[,i])
    d2g_dbeta2_i=crossprod(d2h_dbeta2_i,y/(rates^2))-Matrix(colSums(d2lambda_dbeta2_i))
    H_beta[,i]=d2g_dbeta2_i
    d2h_dudbeta_i=t(t(dh_du)*as.vector(Z[,i]))*as.vector(rates)-dh_du*as.vector(dh_dbeta[,i])
    d2g_dudbeta_i=crossprod(d2h_dudbeta_i,y/(rates^2))-dlambda_dbeta[,i]
    H_ubeta[,i]=d2g_dudbeta_i
  }
  u_ij=lapply(u_nonzeros,h_nonzeros,dh_du,rates)
  h_ij=sparseMatrix(i=unlist(u_nonzeros),
                    j=rep(1:N,sapply(u_nonzeros,length)),
                    x=unlist(u_ij))
  h_ii=crossprod(dh_du,y/rates)-Matrix(colSums(lambda))
  H_u=Diagonal(N,x=as.vector(h_ii))-h_ij
  return( cbind2( rbind2(H_beta,H_ubeta),rbind2(t(H_ubeta),H_u) )-Q_x )
}

eval_fisher<-function(x){
  beta=Matrix(x[1:P],1,P)
  u=Matrix(x[-c(1:P)],1,N)
  lambda=exp(tcrossprod(beta,Z)+u)
  rates=tcrossprod(A,lambda)
  dlambda_dbeta=Z*as.vector(lambda)
  dh_dbeta=A%*%dlambda_dbeta
  dh_du=t(t(A)*as.vector(lambda))
  H_beta=Matrix(0,P,P)
  H_ubeta=Matrix(0,N,P)
  for (i in 1:P){
    d2lambda_dbeta2_i=dlambda_dbeta*as.vector(Z[,i])
    d2h_dbeta2_i=(A%*%d2lambda_dbeta2_i)*as.vector(rates)-dh_dbeta*as.vector(dh_dbeta[,i])
    d2g_dbeta2_i=crossprod(d2h_dbeta2_i,1/rates)-Matrix(colSums(d2lambda_dbeta2_i))
    H_beta[,i]=d2g_dbeta2_i
    d2h_dudbeta_i=t(t(dh_du)*as.vector(Z[,i]))*as.vector(rates)-dh_du*as.vector(dh_dbeta[,i])
    d2g_dudbeta_i=crossprod(d2h_dudbeta_i,1/rates)-dlambda_dbeta[,i]
    H_ubeta[,i]=d2g_dudbeta_i
  }
  u_ij=lapply(u_nonzeros,h_nonzeros,dh_du,rates)
  h_ij=sparseMatrix(i=unlist(u_nonzeros),
                    j=rep(1:N,sapply(u_nonzeros,length)),
                    x=unlist(u_ij))
  h_ii=colSums(dh_du)-colSums(lambda)
  H_u=Diagonal(N,x=h_ii)-h_ij
  return( cbind2( rbind2(H_beta,H_ubeta),rbind2(t(H_ubeta),H_u) ) )
}

is.PD<-function(M){
  tryCatch(Cholesky(M,LDL=F), error=function(e) {NA}, warning=function(w) {NA})
}
near.PD<-function(M,s=1e-04){
  while( isTRUE(class(is.PD(M))=="logical") ){
    s=2*s
    M=M+s*Diagonal(ncol(M))
  }
  return(M)
}

mode_x_fisher<-function(x,k){
  Q_x=bdiag(Q_beta,k*W)
  x0=x
  H0=-eval_fisher(x0)+Q_x
  nablafx0=eval_gradientf(x0,k)
  print(paste("Estimate is",round(x0[1],2),"with",round(max(abs(nablafx0)),2)))
  while( max(abs(nablafx0))>1e-1 ){
    L0=Cholesky(near.PD(H0))
    Delta=solve(L0,nablafx0,system="A")
    xt=x0+Delta
    nablafxt=eval_gradientf(xt,k)
    print(paste("Estimate is",round(xt[1],2),"with",round(max(abs(nablafxt)),2)))
    while( !is.finite(eval_f(xt,k)) | eval_f(xt,k)-eval_f(x0,k)<0 ){
      Delta=0.5*Delta
      xt=x0+Delta
      nablafxt=eval_gradientf(xt,k)
      print(paste("Estimate is",round(xt[1],2),"with",round(max(abs(nablafxt)),2),"and",sqrt(sum(Delta^2))))
      if(max(abs(nablafxt))<1e-1) break
    }
    x0=xt
    H0=-eval_fisher(x0)+Q_x
    nablafx0=eval_gradientf(x0,k)
  }
  return(x0)
}

mode_x_la<-function(k){
  x0=mode_x_fisher(x=mu_x,k=k)
  #x0=mu_x
  H0=-eval_hessianf(x0,k)
  nablafx0=eval_gradientf(x0,k)
  print(paste("Estimate for",round(k,10),"is",round(x0[1],2),"with",round(max(abs(nablafx0)),2)))
  while( max(abs(nablafx0))>1e-1 | isTRUE(class(is.PD(H0))=="logical") ){
    L0=Cholesky(near.PD(H0))
    Delta=solve(L0,nablafx0,system="A")
    xt=x0+Delta
    nablafxt=eval_gradientf(xt,k)
    print(paste("Estimate for",round(k,10),"is",round(xt[1],2),"with",round(max(abs(nablafxt)),2)))
    while( !is.finite(eval_f(xt,k)) | eval_f(xt,k)-eval_f(x0,k)<0 ){
      Delta=0.5*Delta
      xt=x0+Delta
      nablafxt=eval_gradientf(xt,k)
      print(paste("Estimate for",round(k,10),"is",round(xt[1],2),"with",round(max(abs(nablafxt)),2),"and",sqrt(sum(Delta^2))))
      if(max(abs(nablafxt))<1e-1) break
      #if(max(abs(nablafxt))<1e-1 & isFALSE(class(is.PD(-eval_hessianf(xt,k)))=="logical")) break
    }
    x0=xt
    H0=-eval_hessianf(x0,k)
    nablafx0=eval_gradientf(x0,k)
  }
  return(x0)
}

### Approximation of the posterior marginal of k
## Logarithm of the marginal posterior of k
# For mu_mle=c(10,rep(0,N)), mu_beta=10 and Q_beta=1e-2: 
approx_marginalk_la<-function(k){
  x_star=mode_x_la(k=k)
  Q=forceSymmetric(-eval_hessianf(x_star,k))
  L=Cholesky(Q,LDL=FALSE)
  v=eval_gradientf(x_star,k)+crossprod(Q,x_star)
  u=solve(L,v,system="A")
  logn=dgamma(k,shape=a_k+0.5*(N-1),rate=b_k+0.5*q_u,log=T)
  logd=as.numeric(determinant(L)$modulus)-0.5*as.numeric(crossprod(v,u))
  return(logn-logd-457788900)
}


tic<-Sys.time()
approx_marginalk_la(k=1e-5)
toc<-Sys.time()


library(doParallel)
cl<-makeCluster(7)
registerDoParallel(cl)
tic<-Sys.time()
optimize_k<-foreach(k=seq(7.1149e-4,7.1155e-4,by=1e-8),.combine='cbind',.packages=c('Matrix')) %dopar% {
  logmarginal_k=approx_marginalk_la(k)
  return( logmarginal_k )
}
stopCluster(cl)
toc<-Sys.time()



k_star<-7.1152e-4
x_star<-mode_x_la(k=k_star)
save(x_star,file="x_star.RData")
load("./x_star.RData")
H_star<-(-eval_hessianf(x=x_star,k=k_star))
L_star<-Cholesky(H_star,LDL=F)

## x samples
S<-1000
set.seed(1987)
z<-Matrix(rnorm((N+P)*S),nrow=N+P,ncol=S)
s_samples<-solve(L_star,z,system="Lt")
x_samples<-solve(L_star,s_samples,system="Pt")+x_star
T_lambda<-cbind2(Z,Diagonal(N))
loglambda_samples<-T_lambda%*%x_samples
save(loglambda_samples,file="loglambda_samples.RData")


Areas<-st_intersection(BE_LAU2,BE_grid) %>%
  mutate(GC_ID=as.character(GC_ID),LAU2_ID=as.character(LAU2_ID)) %>%
  mutate(Area_km2=as.numeric(st_area(geometry))/1000000) %>%
  st_set_geometry(NULL) %>%
  mutate(GC_ID=as.factor(GC_ID),LAU2_ID=as.factor(LAU2_ID)) %>%
  mutate(GC=as.numeric(GC_ID),LAU2=as.numeric(LAU2_ID)) %>%
  select(GC,LAU2,Area_km2)



logintensity<-function(region,loglambda){
  gridcells=which(A[region,]!=0)
  lx=loglambda[gridcells]+log(A[region,gridcells])
  return(logSumExp(lx))
}

library(doParallel)

cl<-makeCluster(10)
registerDoParallel(cl)
logrates<-foreach(s=1:S,.combine='cbind',.packages=c('Matrix','matrixStats')) %dopar% {
  logrates=apply(array(1:M),1,logintensity,loglambda=loglambda_samples[,s])
    return(logrates)
}
stopCluster(cl)

library(doRNG)
cl<-makeCluster(10)
registerDoRNG(2021)
registerDoParallel(cl)
ymn_samples<-foreach(i=1:nrow(Areas),.combine='rbind',.packages=c('Matrix')) %dopar% {
  samples<-rbinom(S,
                  size=y[Areas$LAU2[i]],
                  prob=exp(log(Areas$Area_km2[i])+loglambda_samples[Areas$GC[i],]-logrates[Areas$LAU2[i],]))
  return(c(GC=Areas$GC[i],LAU2=Areas$LAU2[i],samples))
  
}
stopCluster(cl)

latent_samples<-left_join(Areas,ymn_samples,by=c("GC","LAU2"),copy=TRUE) %>%
  group_by(GC) %>%
  summarise(across(.cols=starts_with("V"),.fns=sum)) %>% 
  select(starts_with("V"))
save(latent_samples,file="latent_samples.RData")
