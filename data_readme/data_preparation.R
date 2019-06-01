setwd("data_readme/raw_data")
library(sf)
library(dplyr)

#Shape file of edges (a line per edge)
rds<-st_read("TRONCON_ROUTE.shp")
#Create a speed column according to road attributes
rds$vitesse<-NA
rds[rds$VOCATION=="Type autoroutier" & rds$ACCES=="A péage","vitesse"]<-110
rds[rds$VOCATION=="Type autoroutier" & rds$ACCES=="Libre","vitesse"]<-90
rds[rds$VOCATION=="Liaison principale","vitesse"]<-70
rds[rds$VOCATION=="Liaison régionale","vitesse"]<-60
rds[rds$VOCATION=="Liaison locale","vitesse"]<-50
#Transformation distance -> travel time
rds$cost<-60*rds$LONGUEUR/rds$vitesse

#Extraction of the starting point of each edge
pts<-st_line_sample(rds,sample=0)
#Extraction of the ending point of each edge
pts2<-st_line_sample(rds,sample=1)
#Merge
ptstot<-(c(pts,pts2))

#Get coordinates of all nodes
ptstotcoord<-st_coordinates(ptstot)
#Remove duplicates
ptstotcoord<-ptstotcoord[!duplicated(ptstotcoord[,1:2]),]
#Name each node with an unique ID
ptstotcoord<-cbind(ptstotcoord,ID=0:(nrow(ptstotcoord)-1))
#Get coordinates of starting nodes
ptscoord<-st_coordinates(pts)
#Get coordinates of ending nodes
pts2coord<-st_coordinates(pts2)
ptscoord<-data.frame(ptscoord)
#matrix to data frame
pts2coord<-data.frame(pts2coord)
ptstotcoord<-data.frame(ptstotcoord)
#Merge unique ID to starting nodes
ptscoord<-left_join(ptscoord,ptstotcoord[,c(1,2,5)],by=c("X"="X","Y"="Y"))
#Merge unique ID to ending nodes
pts2coord<-left_join(pts2coord,ptstotcoord[,c(1,2,5)],by=c("X"="X","Y"="Y"))
colnames(ptscoord)[5]<-"point1"
colnames(pts2coord)[5]<-"point2"
#Add nodes to the edges data
rds<-cbind(rds,point1=ptscoord$point1,point2=pts2coord$point2)

##Create final symbolic graph : from,to, weight
#Undirected eges and directed edges
options(scipen=999)
dat<-data.frame(depart=rds$point1,arrivee=rds$point2,cost=rds$cost)
dat<-dat[rds$SENS %in% c("Sens unique","Double sens"),]
#Undirected eges and opposed directed edges
dat2<-data.frame(depart=rds$point2,arrivee=rds$point1,cost=rds$cost)
dat2<-dat2[rds$SENS %in% c("Sens inverse","Double sens"),]
dat<-rbind(dat,dat2)
dat<-dat[dat$depart!=dat$arrivee,]
colnames(dat)<-c("from","to","dist")
write.csv(dat,"data_readme/roads.csv",row.names = F)