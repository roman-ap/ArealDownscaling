##### Discrepancy measures
setwd("~/Nextcloud/UofE/Project/Belgium")


load("./BE_LAU2.RData")
load("./BE_LAU2pop.RData")
load("./BE_grid.RData")
load("./BE_gridpop.RData")
load("./k_weights1.RData")
load("./k_weights2.RData")
load("./k_weights3.RData")
load("./k_weights4.RData")
load("./x_means1.RData")
load("./x_means2.RData")
load("./x_means3.RData")
load("./x_means4.RData")
load("./x_vars1.RData")
load("./x_vars2.RData")
load("./x_vars3.RData")
load("./x_vars4.RData")
load("./loglambda_samples.RData")
load("./latent_samples.RData")


library(sf)
library(dplyr)
library(ggplot2)
library(latex2exp)


N<-31512L
P<-2
means1<-sapply(1:N,function(n) weighted.mean(x_means1[n,],w=k_weights1))
means2<-sapply(1:N,function(n) weighted.mean(x_means2[n,],w=k_weights2))
means3<-sapply(1:(N+P),function(n) weighted.mean(x_means3[n,],w=k_weights3))
means4<-sapply(1:(N+P),function(n) weighted.mean(x_means4[n,],w=k_weights4))
means_poi<-latent_samples %>% mutate(Means=rowMeans(.)) %>% select(Means)
vars1<-sapply(1:N,function(n) weighted.mean((x_means1[n,])^2 + x_vars1[n,] - means1[n]^2,w=k_weights1))
vars2<-sapply(1:N,function(n) weighted.mean((x_means2[n,])^2 + x_vars2[n,] - means2[n]^2,w=k_weights2))
vars3<-sapply(1:(N+P),function(n) weighted.mean((x_means3[n,])^2 + x_vars3[n,] - means3[n]^2,w=k_weights3))
vars4<-sapply(1:(N+P),function(n) weighted.mean((x_means4[n,])^2 + x_vars4[n,] - means4[n]^2,w=k_weights4))



BE_aw<-left_join(BE_LAU2,BE_LAU2pop,by="LAU2_ID") %>% 
  mutate(Area_km2=as.numeric(st_area(geometry))/1000000,aw_est=Population/Area_km2) %>%
  select(LAU2_ID,aw_est) %>%
  st_set_geometry(NULL)
BEgridded_LAU2<-st_intersection(BE_grid,BE_LAU2) %>% 
  mutate(Area_km2=as.numeric(st_area(GEOMETRY))/1000000) %>%
  st_set_geometry(NULL)
BEgridded_aw<-left_join(BEgridded_LAU2,BE_aw,by="LAU2_ID") %>% 
  mutate(est=Area_km2*aw_est) %>%
  group_by(GC_ID) %>%
  summarise(Pop_est=sum(est)) %>%
  left_join(BE_gridpop,by="GC_ID") %>%
  mutate(Pop=if_else(is.na(TOT_P),0L,TOT_P),Delta=abs(Pop_est-Pop))
Delta_aw<-0.5*sum(BEgridded_aw$Delta)


BEgridded_gmrf1<-left_join(BE_grid,BE_gridpop,by="GC_ID") %>% 
  mutate(Pop=if_else(is.na(TOT_P),0L,TOT_P),Pop_est=means1,Delta=abs(Pop_est-Pop),
         p_values=pnorm(Pop,mean=means1,sd=sqrt(vars1)))

Delta_gmrf1<-0.5*sum(BEgridded_gmrf1$Delta)

ggplot(BEgridded_gmrf1,aes(x=p_values))+
  geom_histogram(binwidth=0.01,color="black",fill="white")+
  expand_limits(x=c(0,1))+
  labs(x=TeX("$p_{i}$"))+
  theme_bw()


BEgridded_gmrf2<-left_join(BE_grid,BE_gridpop,by="GC_ID") %>% 
  mutate(Pop=if_else(is.na(TOT_P),0L,TOT_P),Pop_est=means2,Delta=abs(Pop_est-Pop),
         p_values=pnorm(Pop,mean=means2,sd=sqrt(vars2)))

Delta_gmrf2<-0.5*sum(BEgridded_gmrf2$Delta) 

ggplot(BEgridded_gmrf2,aes(x=p_values))+
  geom_histogram(binwidth=0.01,color="black",fill="white")+
  expand_limits(x=c(0,1))+
  labs(x=TeX("$p_{i}$"))+
  theme_bw()

BEgridded_gmrf3<-left_join(BE_grid,BE_gridpop,by="GC_ID") %>% 
  mutate(Pop=if_else(is.na(TOT_P),0L,TOT_P),Pop_est=as.vector(Tx%*%means3),Delta=abs(Pop_est-Pop),
         p_values=pnorm(Pop,mean=as.vector(Tx%*%means3),sd=as.vector(sqrt(Tx%*%vars3))))

Delta_gmrf3<-0.5*sum(BEgridded_gmrf3$Delta)

ggplot(BEgridded_gmrf3,aes(x=p_values))+
  geom_histogram(binwidth=0.01,color="black",fill="white")+
  expand_limits(x=c(0,1))+
  labs(x=TeX("$p_{i}$"))+
  theme_bw()

BEgridded_gmrf4<-left_join(BE_grid,BE_gridpop,by="GC_ID") %>% 
  mutate(Pop=if_else(is.na(TOT_P),0L,TOT_P),Pop_est=as.vector(Tx%*%means4),Delta=abs(Pop_est-Pop),
         p_values=pnorm(Pop,mean=as.vector(Tx%*%means4),sd=as.vector(sqrt(Tx%*%vars4))))

Delta_gmrf4<-0.5*sum(BEgridded_gmrf4$Delta) 

ggplot(BEgridded_gmrf4,aes(x=p_values))+
  geom_histogram(binwidth=0.01,color="black",fill="white")+
  expand_limits(x=c(0,1))+
  labs(x=TeX("$p_{i}$"))+
  theme_bw()


BEgridded_poi<-left_join(BE_grid,BE_gridpop,by="GC_ID") %>% 
  mutate(Pop=if_else(is.na(TOT_P),0L,TOT_P),Pop_est=means_poi,Delta=abs(Pop_est-Pop) )

Delta_poi<-0.5*sum(BEgridded_poi$Delta)

maximo<-data.frame(muestra=as.matrix(latent_samples)[11478,])

jpeg(file="/home/koko/Nextcloud/UofE/Documents/Thesis/Figures/sample_11478.jpeg",width=600,height=300)
ggplot(maximo,aes(x=muestra))+
  geom_histogram(bins=100,color="black",fill="white")+
  geom_vline(xintercept=BEgridded_poi$Pop[11478],color="red")+
  labs(x=TeX("$\\widetilde{y}_{11478}$"))+
  theme_bw()
dev.off()
