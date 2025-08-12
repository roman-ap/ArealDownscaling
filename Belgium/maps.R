########## Plots
setwd("~/Nextcloud/UofE/Project/Belgium")

### Datasets
load("./BE_LAU2.RData")
load("./BE_LAU2pop.RData")
load("./BE_grid.RData")
load("./BE_gridpop.RData")
load("./BE_CLC.RData")
load("./BE.RData")
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



### Packages
library(sf)
library(dplyr)
library(latex2exp)
library(ggplot2)
library(ggspatial)
library(viridis)


N<-31512L
P<-2
means1<-sapply(1:N,function(n) weighted.mean(x_means1[n,],w=k_weights1))
means2<-sapply(1:N,function(n) weighted.mean(x_means2[n,],w=k_weights2))
means3<-sapply(1:(N+P),function(n) weighted.mean(x_means3[n,],w=k_weights3))
means4<-sapply(1:(N+P),function(n) weighted.mean(x_means4[n,],w=k_weights4))
loglambda_means<-apply(loglambda_samples,1,mean)
latent_means<-apply(latent_samples,1,mean)
vars1<-sapply(1:N,function(n) weighted.mean((x_means1[n,])^2 + x_vars1[n,] - means1[n]^2,w=k_weights1))
vars2<-sapply(1:N,function(n) weighted.mean((x_means2[n,])^2 + x_vars2[n,] - means2[n]^2,w=k_weights2))
vars3<-sapply(1:(N+P),function(n) weighted.mean((x_means3[n,])^2 + x_vars3[n,] - means3[n]^2,w=k_weights3))
vars4<-sapply(1:(N+P),function(n) weighted.mean((x_means4[n,])^2 + x_vars4[n,] - means4[n]^2,w=k_weights4))
loglambda_vars<-apply(loglambda_samples,1,var)
latent_vars<-apply(latent_samples,1,var)


BElau2<-left_join(BE_LAU2,BE_LAU2pop,by='LAU2_ID') %>%
  mutate(Area=as.numeric(st_area(geometry)/1000000),Density=Population/Area)

BEgrid<-left_join(BE_grid,BE_gridpop,by='GC_ID') %>%
  rename(True=TOT_P) %>%
  mutate(#M1=means1,
         #M2=means2,
         #M3=as.vector(Tx%*%means3),
         #M4=as.vector(Tx%*%means4),
         LM=loglambda_means,
         M5=latent_means,
         #SD1=sqrt(vars1),
         #SD2=sqrt(vars2),
         #SD3=as.vector(sqrt(Tx%*%vars3)),
         #SD4=as.vector(sqrt(Tx%*%vars4)),
         LSD=sqrt(loglambda_vars),
         SD5=sqrt(latent_vars)
         )
  
BE211lau2<-filter(BElau2,stringr::str_sub(LAU2_ID,1,5)=="BE211")
BE211<-st_union(BE211lau2)
BE211_match<-st_intersects(BE211,BEgrid)
BE211grid<-BEgrid[unlist(BE211_match),]



ggplot(data=BE211grid) + 
  geom_errorbar(mapping=aes(x=GC_ID, ymin=Q05Est1, ymax=Q95Est1),width=0.2,colour="blue") + 
  geom_point(mapping=aes(x=GC_ID,y=MeanEst1),shape=21, fill="white") +
  geom_point(mapping=aes(x=GC_ID,y=True),shape=21, fill="red") +
  theme(panel.background=element_blank(),legend.title=element_blank())



ggplot()+
  geom_sf(data=BE211lau2,color="black",aes(fill=Density))+
  scale_fill_gradient(low="green",high="red",limits=c(0,4000),name="Inhabitants")+
  annotation_scale(location="br",width_hint=0.1)+
  annotation_north_arrow(location="tr",which_north="true",height=unit(1,"cm"),width=unit(1,"cm"))+
  coord_sf(crs=st_crs(3035),datum=NA)+
  theme(panel.background=element_blank(),legend.position="right")

ggplot()+
  geom_sf(data=BE211grid,aes(fill=True),color='gray25',size=0.1)+
  scale_fill_gradient(low="green",high="red",limits=c(0,20000),na.value="white",name="Inhabitants")+
  geom_sf(data=BE211lau2,color="black",fill=NA)+
  annotation_scale(location="br",width_hint=0.1)+
  annotation_north_arrow(location="tr",which_north="true",height=unit(1,"cm"),width=unit(1,"cm"))+
  coord_sf(crs=st_crs(3035),datum=NA)+
  theme(panel.background=element_blank(),legend.position="right")



ggplot()+
  geom_sf(data=BE211grid,aes(fill=M1),color='gray25',size=0.1)+
  scale_fill_gradientn(colours=c("white","white","green","red"),
                       values=c(0,0.002/4130.001,0.0021/4130.001,1),
                       limits=c(-0.001,4130),
                       na.value="darkgray",
                       name="Inhabitants")+
  geom_sf(data=BE211lau2,color="black",fill=NA)+
  annotation_scale(location="br",width_hint=0.1)+
  annotation_north_arrow(location="tr",which_north="true",height=unit(1,"cm"),width=unit(1,"cm"))+
  coord_sf(crs=st_crs(3035),datum=NA)+
  theme(panel.background=element_blank(),legend.position="right")

ggplot()+
  geom_sf(data=BE211grid,aes(fill=SD1),color='gray25',size=0.1)+
  scale_fill_viridis(direction=-1,option="magma")+
  geom_sf(data=BE211lau2,color="black",fill=NA)+
  annotation_scale(location="br",width_hint=0.1)+
  annotation_north_arrow(location="tr",which_north="true",height=unit(1,"cm"),width=unit(1,"cm"))+
  coord_sf(crs=st_crs(3035),datum=NA)+
  theme(panel.background=element_blank(),legend.position="right",legend.title=element_blank())



ggplot()+
  geom_sf(data=BE211grid,aes(fill=M2),color='gray25',size=0.1)+
  scale_fill_gradientn(colours=c("white","white","green","red"),
                       values=c(0,0.2/8700.1,0.21/8700.1,1),
                       limits=c(-0.1,8700),
                       na.value="darkgray",
                       name="Inhabitants")+
  geom_sf(data=BE211lau2,color="black",fill=NA)+
  annotation_scale(location="br",width_hint=0.1)+
  annotation_north_arrow(location="tr",which_north="true",height=unit(1,"cm"),width=unit(1,"cm"))+
  coord_sf(crs=st_crs(3035),datum=NA)+
  theme(panel.background=element_blank(),legend.position="right")

ggplot()+
  geom_sf(data=BE211grid,aes(fill=SD2),color='gray25',size=0.1)+
  scale_fill_viridis(direction=-1,option="magma")+
  geom_sf(data=BE211lau2,color="black",fill=NA)+
  annotation_scale(location="br",width_hint=0.1)+
  annotation_north_arrow(location="tr",which_north="true",height=unit(1,"cm"),width=unit(1,"cm"))+
  coord_sf(crs=st_crs(3035),datum=NA)+
  theme(panel.background=element_blank(),legend.position="right",legend.title=element_blank())



ggplot()+
  geom_sf(data=BE211grid,aes(fill=M3),color='gray25',size=0.1)+
  scale_fill_gradientn(colours=c("white","white","green","red"),
                       values=c(0,0.2/20300.1,0.21/20300.1,1),
                       limits=c(-0.1,20300),
                       na.value="darkgray",
                       name="Inhabitants")+
  geom_sf(data=BE211lau2,color="black",fill=NA)+
  annotation_scale(location="br",width_hint=0.1)+
  annotation_north_arrow(location="tr",which_north="true",height=unit(1,"cm"),width=unit(1,"cm"))+
  coord_sf(crs=st_crs(3035),datum=NA)+
  theme(panel.background=element_blank(),legend.position="right")

ggplot()+
  geom_sf(data=BE211grid,aes(fill=SD3),color='gray25',size=0.1)+
  scale_fill_viridis(direction=-1,option="magma")+
  annotation_scale(location="br",width_hint=0.1)+
  annotation_north_arrow(location="tr",which_north="true",height=unit(1,"cm"),width=unit(1,"cm"))+
  geom_sf(data=BE211lau2,color="black",fill=NA)+
  coord_sf(crs=st_crs(3035),datum=NA)+
  theme(panel.background=element_blank(),legend.position="right",legend.title=element_blank())



ggplot()+
  geom_sf(data=BE211grid,aes(fill=M4),color='gray25',size=0.1)+
  scale_fill_gradientn(colours=c("white","white","green","red"),
                       values=c(0,0.2/21400.1,0.21/21400.1,1),
                       limits=c(-0.1,21400),
                       na.value="darkgray",
                       name="Inhabitants")+
  geom_sf(data=BE211lau2,color="black",fill=NA)+
  annotation_scale(location="br",width_hint=0.1)+
  annotation_north_arrow(location="tr",which_north="true",height=unit(1,"cm"),width=unit(1,"cm"))+
  coord_sf(crs=st_crs(3035),datum=NA)+
  theme(panel.background=element_blank(),legend.position="right")

ggplot()+
  geom_sf(data=BE211grid,aes(fill=SD4),color='gray25',size=0.1)+
  scale_fill_viridis(direction=-1,option="magma")+
  annotation_scale(location="br",width_hint=0.1)+
  annotation_north_arrow(location="tr",which_north="true",height=unit(1,"cm"),width=unit(1,"cm"))+
  geom_sf(data=BE211lau2,color="black",fill=NA)+
  coord_sf(crs=st_crs(3035),datum=NA)+
  theme(panel.background=element_blank(),legend.position="right",legend.title=element_blank())


jpeg(file="/home/koko/Nextcloud/UofE/Documents/Thesis/Figures/BE211_loglambda_mean.jpeg",width=900,height=900)
ggplot()+
  geom_sf(data=BE211grid,aes(fill=LM),color='gray25',size=0.1)+
  scale_fill_gradient(low="grey",high="brown",name=TeX("$E( \\log \\lambda_{n} )$"))+
  geom_sf(data=BE211lau2,color="black",fill=NA)+
  annotation_scale(location="br",width_hint=0.1)+
  annotation_north_arrow(location="tr",which_north="true",height=unit(1,"cm"),width=unit(1,"cm"))+
  coord_sf(crs=st_crs(3035),datum=NA)+
  theme(panel.background=element_blank(),legend.position="right")
dev.off()

jpeg(file="/home/koko/Nextcloud/UofE/Documents/Thesis/Figures/BE211_loglambda_sd.jpeg",width=900,height=900)
ggplot()+
  geom_sf(data=BE211grid,aes(fill=LSD),color='gray25',size=0.1)+
  scale_fill_viridis(direction=-1,option="magma",name=TeX("$\\sigma( \\log \\lambda_{n} )$"))+
  annotation_scale(location="br",width_hint=0.1)+
  annotation_north_arrow(location="tr",which_north="true",height=unit(1,"cm"),width=unit(1,"cm"))+
  geom_sf(data=BE211lau2,color="black",fill=NA)+
  coord_sf(crs=st_crs(3035),datum=NA)+
  theme(panel.background=element_blank(),legend.position="right",legend.title=element_blank())
dev.off()


ggplot()+
  geom_sf(data=BE211grid,aes(fill=M5),color='gray25',size=0.1)+
  scale_fill_gradient(low="green",high="red",limits=c(1,10850),na.value="white",name="Inhabitants")+
  geom_sf(data=BE211lau2,color="black",fill=NA)+
  annotation_scale(location="br",width_hint=0.1)+
  annotation_north_arrow(location="tr",which_north="true",height=unit(1,"cm"),width=unit(1,"cm"))+
  coord_sf(crs=st_crs(3035),datum=NA)+
  theme(panel.background=element_blank(),legend.position="right")

ggplot()+
  geom_sf(data=BE211grid,aes(fill=SD5),color='gray25',size=0.1)+
  scale_fill_viridis(direction=-1,option="magma")+
  annotation_scale(location="br",width_hint=0.1)+
  annotation_north_arrow(location="tr",which_north="true",height=unit(1,"cm"),width=unit(1,"cm"))+
  geom_sf(data=BE211lau2,color="black",fill=NA)+
  coord_sf(crs=st_crs(3035),datum=NA)+
  theme(panel.background=element_blank(),legend.position="right",legend.title=element_blank())

  
  
  
BE211clc_match<-st_intersects(BE211,BE_CLC)
BE211_CLC<-BE_CLC[unlist(BE211clc_match),] %>%
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
    st_intersection(BE211)
  
ggplot()+
  geom_sf(data=BE211_CLC,aes(fill=myCLC_class),color='gray25',size=0.1)+
  scale_fill_manual(values=c("Urban"="darkred",
                             "Rural"="orangered",
                             "Commercial"="grey75",
                             "Transport"="grey50",
                             "Sites"="grey25",
                             "Leisure"="pink",
                             "Arable"="gold",
                             "Crops"="wheat",
                             "Pasture"="green",
                             "Groves"="peru",
                             "Forest"="forestgreen",
                             "Vegetation"="darkolivegreen",
                             "Bare"="khaki",
                             "Wetlands"="tan",
                             "Water"="dodgerblue"),na.value="white",name="Category")+
  geom_sf(data=BE211lau2,color="black",fill=NA)+
  annotation_scale(location="br",width_hint=0.1)+
  annotation_north_arrow(location="tr",which_north="true",height=unit(1,"cm"),width=unit(1,"cm"))+
  coord_sf(crs=st_crs(3035),datum=NA)+
  theme(panel.background=element_blank(),legend.position="right")



BE21111056<-filter(BElau2,LAU2_ID=="BE211_11056")
BE21111056grid_match<-st_intersects(BE21111056,BEgrid)
BE21111056grid<-BEgrid[unlist(BE21111056grid_match),]
BE21111056clc_match<-st_intersects(BE21111056grid,BE_CLC)
BE21111056clc<-BE_CLC[unlist(BE21111056clc_match),] %>%
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
  st_intersection(BE21111056grid)
myBE21111056<-filter(BE,GC_ID %in% BE21111056grid$GC_ID) %>%
  st_set_geometry(NULL) %>%
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
  left_join(BE21111056grid,by="GC_ID") %>%
  select(GC_ID,True,Rural)

library(Matrix)
library(igraph)
BE21111056centroids<-st_centroid(st_geometry(BE21111056grid))
BE21111056neigbours<-1*st_relate(BE21111056grid,BE21111056grid,pattern="F***1****",sparse=FALSE)
BE21111056graph<-graph_from_adjacency_matrix(BE21111056neigbours,mode="upper")
zeros<-if_else(myBE21111056$Rural>0,F,T)
zeros_color<-if_else(zeros,"red","black")
grid_match<-st_relate(BE21111056grid,BE21111056grid,pattern="F***1****",sparse=TRUE)
nonzeros_match<-grid_match[1:37]
nonzeros_match[which(zeros==TRUE)]<-lapply(grid_match[zeros==TRUE],setdiff,y=which(zeros==TRUE)) 
nonzeros_neigbours<-sparseMatrix(i=c( rep(1:37,sapply(nonzeros_match,length)) ),
                                 j=c( unlist(nonzeros_match) ),
                                 x=c( rep(1,length(unlist(nonzeros_match))) ),
                                 dims=c(37,37))
BE21111056graph_nonzeros<-graph_from_adjacency_matrix(nonzeros_neigbours,mode="upper")

ggplot()+
  geom_sf(data=BE21111056clc,aes(fill=myCLC_class),color='gray25',size=0.1)+
  scale_fill_manual(values=c("Urban"="darkred",
                             "Rural"="orangered",
                             "Commercial"="grey75",
                             "Transport"="grey50",
                             "Sites"="grey25",
                             "Leisure"="pink",
                             "Arable"="gold",
                             "Crops"="wheat",
                             "Pasture"="green",
                             "Groves"="peru",
                             "Forest"="forestgreen",
                             "Vegetation"="darkolivegreen",
                             "Bare"="khaki",
                             "Wetlands"="tan",
                             "Water"="dodgerblue"),na.value="white",name="Category")+
  geom_sf(data=BE21111056,color="black",lwd=1,fill=NA)+
  geom_sf(data=BE21111056grid,color="black",lwd=0.1,fill=NA)+
  annotation_scale(location="br",width_hint=0.25)+
  annotation_north_arrow(location="tr",which_north="true",height=unit(0.75,"cm"),width=unit(0.75,"cm"))+
  coord_sf(crs=st_crs(3035),datum=NA)+
  theme(panel.background=element_blank(),legend.position="right")
par(mfrow=c(1,2),mar=c(0,0,0,0))
plot(BE21111056graph,vertex.size=5,vertex.color=zeros_color,vertex.label=NA,edge.color="black",layout=st_coordinates(BE21111056centroids))
plot(BE21111056graph_nonzeros,vertex.size=5,vertex.color=zeros_color,vertex.label=NA,edge.color="black",layout=st_coordinates(BE21111056centroids))


