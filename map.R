#install.packages("rworldmap")
library(rworldmap)
library(readr)
library(stringr)
library(RgoogleMaps)
#MDSCoord_1 <- read_csv("Coordonnées_MDS_dim123/MDSCoord_1.csv")


####### Clustering de l clusters
fit2 = cmdscale(jaccard_abundance ,eig=TRUE, k=10)
#gps=GPScoordinates2 #"Lat"  "Long"
gps=design[,c(1,4,5)]
names(gps)
l = 2
res2 = kmeans(fit2$points,l, nstart = 1000)$cluster


# map("worldHires")
# clus_colors<-rainbow(8)
# points(gps$Long,gps$Lat,pch=16,cex=1)

newmap <- getMap(resolution = "li")

color = c()
labels = c()
cl = res2 # clustering

for(i in 1:length(cl)){
  if(cl[i]==1){
    color[i] = "blue"
  }
  if(cl[i]==2){
    color[i] = "red"
  }
  if(cl[i]==3){
    color[i] = "gold4"
  }
  if(str_length(design$Station[i]) == 5){
    labels[i] = substr(design$Station[i],1,1)
  }
  if(str_length(design$Station[i]) == 6){
    labels[i] = substr(design$Station[i],1,2)
  }
  if(str_length(design$Station[i]) == 7){
    labels[i] = substr(design$Station[i],1,3)
  }
}

DCM_indices = which(str_detect(design$Station,"DCM")==TRUE)
SUR_indices = which(str_detect(design$Station,"SUR")==TRUE)

SUR_only_indices = c()
a = 1
for(i in 1:length(cl)){
  if(length(which(labels[SUR_indices][i]!=labels[DCM_indices]))==length(labels[DCM_indices])){
    SUR_only_indices[a] = SUR_indices[i]
    a = a + 1
  }
}

plot(newmap, xlim = c(-180, 90), ylim = c(-75, 75), asp = 1)
points(gps$Long[DCM_indices], gps$Lat[DCM_indices], col = color[DCM_indices], cex = 1, pch = 17)
points(gps$Long[SUR_indices], gps$Lat[SUR_indices] + 2, col = color[SUR_indices], cex = 1, pch = 16)
text(gps$Long[DCM_indices], gps$Lat[DCM_indices], labels = labels[DCM_indices], pos = 1, col = "darkgreen", cex = 0.8, lty = 1)

text(gps$Long[SUR_only_indices], gps$Lat[SUR_only_indices] + 2, labels = labels[SUR_only_indices], pos = 1, col = "darkgreen", cex = 0.8, lty = 1)
