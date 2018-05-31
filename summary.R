library(vegan)
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(gdata)

#------------------------
data.wd <- "/home/diabakhate/HydroGen_Abdou"
setwd(dir=data.wd) 
#--
data1=Genocenoses_env_parameters_all_tara
names(data1)
dim(data1)
str(data1)
unique(data1$Station)   # On a 159 stations
unique(data1$Fraction)

#---------------------Extraction des différents jeu de données pour chaque fraction de taille
data.0_0.2 = data1[which((as.character(data1$Fraction) == "0-0.2") == TRUE),]    #dim(70,20)
data.0.22_3 = data1[which((as.character(data1$Fraction) == "0.22-3") == TRUE),]  #dim(105,20)
data.0.8_5 = data1[which((as.character(data1$Fraction) == "0.8-5") == TRUE),]    #dim(131,20)
data.5_20 = data1[which((as.character(data1$Fraction) == "43952") == TRUE),]     #dim(77,20)
data.20_180 = data1[which((as.character(data1$Fraction) == "20-180") == TRUE),]  #dim(126,20)
data.180_2000 = data1[which((as.character(data1$Fraction) == "180-2000") == TRUE),] #dim(135,20)

#---------Le nom des stations correspondants à chaque fraction de taille
rownames(data.0_0.2) = data.0_0.2$Station
rownames(data.0.22_3) = data.0.22_3$Station
rownames(data.0.8_5) = data.0.8_5$Station
rownames(data.5_20) = data.5_20$Station
rownames(data.20_180) = data.20_180$Station
rownames(data.180_2000) = data.180_2000$Station

#__________________________________________________________________________________
design = data.0_0.2 # On s'intéresse à la taille de fraction 0_0.2
#_________________________________________________________________________________
names(design)
dim(design)

#--------On s'intéresse aux matrices distances pour la taille de fraction donnée: 0_0.2.....
size_fraction1 = "0-0.2"
size_fraction2 = "0.22-3"
size_fraction3 = "0.8-5"
size_fraction4 = "43952" ########c'est le 5-20
size_fraction5 = "20-180"
size_fraction6 = "180-2000"

#### Compiler la fonction utils d'abord
import_data(size_fraction1, samples = rownames(design))  
# #ATTENTION : bien prendre le bon "size_fraction" associé au bon design (dans ce cas, c'est size_fraction1)

design <- design[metagenomic_sample, ]  #nouveau jeu de taille (70,20)
dim(design)
#rownames(design)==rownames(jaccard_abundance) # pour verifier le bon ordonnancement des noms des stations
names(design)
summary(design)

library(mice)#les NA
mice = complete(mice(design,method="norm.predict",m=1)) # fonctionne dans le cas de valeurs numeriques uniquement

design$T  = mice$T
design$Sal  = mice$Sal
design$Chla  = mice$Chla 
design$O2  = mice$O2
design$NO3   = mice$NO3 
design$Phos   = mice$Phos 
design$Si   = mice$Si 

summary(design)

#-------------

library(MDMR)
D <- dist(jaccard_abundance, method = "euclidean")

new_des=design[,-c(1,2,20)]
names(new_des)
plot(new_des)

design$Genocenose=as.factor(design$Genocenose)
design$month=as.factor(design$month)
str(new_des)

summary(new_des[,-c(1,17)])
new_des_CR=apply(new_des[,-c(1,17)],2,scale)
summary(new_des_CR)
class(new_des_CR)

plot(as.data.frame(new_des_CR))


# design$O2=log10(design$O2)
# design$NO3_m=log10(design$NO3_m)
# design$NO3=log10(design$NO3)
# design$NO2=log10(design$NO2)
# design$NH4=log10(design$NH4)

mdmr.res <- mdmr(X = new_des_CR,D)
summary(mdmr.res)

#----------------

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# MDMR & Adonis
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#----------- Matrice distance
D = jaccard_abundance
n = dim(D)[1]

names(as.data.frame(new_des_CR))
dim(as.data.frame(new_des_CR))
summary(as.data.frame(new_des_CR))

#----X :  CAS D'UN REGRESSEUR 
#parms1 = as.data.frame(new_des_CR)$Long
#X = matrix(parms1,nrow=n,ncol=1)
parms1 = c(design$Lat,design$Long)
X = matrix(parms1,nrow=n,ncol=2)
det(t(X)%*%X) # donc la matrice est inversible

#---H=X.(X'X)^-1.X'
H = X%*%solve(t(X)%*%X)%*%t(X)
H

#----
J = diag(rep(1,n)) - matrix(1,n,n)/n
J
G = -0.5*J%*%D^2%*%J
G

#----------tr(HGH) **************************************************
sum(diag(H%*%G%*%H)) # tr(HGH)=somme(val.propre)  = 32258.88


#----------tr[(I-H)G(I-H)]*****************************************
I=diag(rep(1,n))
Mat=I-H
Mat

sum(diag(Mat%*%G%*%Mat))  # tr[(I-H)G(I-H)]=somme(val.propre)  = 247406.2

#%%%%%%%%%%%%%%%%%% Statistique de test F

adonis(as.dist(jaccard_abundance) ~Lat+Long, data =as.data.frame(new_des_CR))
mdmr.res <- mdmr(X,D)
summary(mdmr.res)

#%%%%%%%%%%%%%%%%%%%%%%%%%% Adonis sur variable qualitative
names(new_des)
new_des$month=as.factor(new_des$month)
#new_des$Genocenose=as.factor(new_des$Genocenose)

adonis(as.dist(jaccard_abundance) ~month, data =new_des)
anova.cca(capscale(as.dist(jaccard_abundance) ~ month, data = new_des),by="terms")
### Pairwise tests
phoc <- with(design, betadisper(as.dist(jaccard_abundance), month))
TukeyHSD(phoc)
plot(phoc)

#with(new_des, anova(capscale(as.dist(jaccard_abundance)~ month, data=new_des), by="term"))
#with(new_des, adonis(as.dist(jaccard_abundance) ~ month, data=new_des, method="euclid"))

#Comparaison avec ADONIS
adonis(as.dist(jaccard_abundance) ~Lat+Long+T+Sal+Chla+O2+NO3_m+NO2+NH4+Fe+SSD+Phos+Si+depth, data =as.data.frame(new_des_CR))
#dbRDA= capscale(jaccard_abundance ~ Lat+Long+T+Sal+Chla+O2+NO3_m+NO2+NH4+Fe+SSD+Phos+Si+depth, data = as.data.frame(new_des_CR), dist="bray")
#anova(dbRDA, by="terms", permu=200)## test for sig.environ. variables

########################################################################

#%%%% Multi Response Permutation Procedure of Within- versus Among-Group Dissimilarities

dune.mrpp <- mrpp(as.dist(jaccard_abundance), new_des$month)
dune.mrpp

#install.packages("BiodiversityR")
#library(BiodiversityR)
#nested.npmanova(as.dist(jaccard_abundance) ~month, method="bray", data =
#                  new_des) 

#------------------------------------------------------------------
library(MASS)
parcoord(new_des_CR, col = as.numeric(new_des$Genocenose), lty = 1)
legend("top", legend=unique(as.numeric(new_des$Genocenose)),
       col=unique(as.numeric(new_des$Genocenose)),lty=1:8,horiz = TRUE, cex=0.8)
title("Parallel Coordinate Plot sur Genocenose")
#------------------------------------------------------------------
#install.packages("MVN")
library(MVN)
result=mvn(new_des_CR,mvnTest = "royston",univariatePlot = "qqplot")
result=mvn(new_des_CR,mvnTest = "royston",univariatePlot = "histogram")
result=mvn(new_des_CR,mvnTest = "royston",univariateTest="SW",desc = TRUE)
result$univariateNormality
result$multivariateNormality
#install.packages("Compositional")






