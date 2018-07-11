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
design = data.0.8_5 # On s'intéresse à la taille de fraction 0_0.2
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
import_data(size_fraction3, samples = rownames(design))  
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
des=new_des[,-c(1,17)]
summary(des)
str(des)


# new_des_CR=apply(new_des[,-c(1,17)],2,scale)
# summary(new_des_CR)
# class(new_des_CR)
# 
# plot(as.data.frame(new_des_CR))


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
str(new_des)
P=adonis(as.dist(jaccard_abundance) ~month, data =new_des)
P
hist(P$f.perms) #To obtain the distribution of pseudo-f's based on H0
###################

# bird.mat=sqrt(jaccard_abundance)
# bird.dist=vegdist(bird.mat,method = 'bray')
# 
# dispersion<-betadisper(bird.dist, group=new_des$month)
# permutest(dispersion)
# plot(dispersion, hull=FALSE, ellipse=TRUE) ##sd ellipse


#################
#anova(betadisper(as.dist(jaccard_abundance), new_des$month),type="centroid")
anova.cca(capscale(as.dist(jaccard_abundance) ~ month, data = new_des),by="terms")
nested.npmanova(jaccard_abundance ~ month, data = new_des)
### Pairwise tests
phoc <- with(new_des, betadisper(as.dist(jaccard_abundance), month))
TukeyHSD(phoc)
plot(phoc)

library(RVAideMemoire)
#pairwise.perm.manova(iris[,1:4],iris$Species,nperm=49)
# or
#pairwise.perm.manova(dist(iris[,1:4],"euclidian"),iris$Species,nperm=49)

#install.packages("PMCMR")
#library(PMCMR)
#kruskal.test(jaccard_abundance,new_des$month)

#with(new_des, anova(capscale(as.dist(jaccard_abundance)~ month, data=new_des), by="term"))
#with(new_des, adonis(as.dist(jaccard_abundance) ~ month, data=new_des, method="euclid"))

#Comparaison avec ADONIS
adonis(as.dist(jaccard_abundance) ~Lat+Long+T+Sal+Chla+O2+NO3_m+NO2+NH4+Fe+SSD+Phos+Si+depth ,data =as.data.frame(new_des_CR),distance="euclidean")
#dbRDA= capscale(jaccard_abundance ~ Lat+Long+T+Sal+Chla+O2+NO3_m+NO2+NH4+Fe+SSD+Phos+Si+depth, data = as.data.frame(new_des_CR), dist="bray")
#anova(dbRDA, by="terms", permu=200)## test for sig.environ. variables

########################################################################

#%%%% Multi Response Permutation Procedure of Within- versus Among-Group Dissimilarities

dune.mrpp <- mrpp(as.dist(jaccard_abundance), new_des$month,distance = "euclidean")
dune.mrpp
#historgram of permuted delta scores
hist(dune.mrpp$boot.deltas)

temp <- meandist(as.dist(jaccard_abundance), new_des$month)
temp
plot(temp)
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
names(new_des)
result=mvn(new_des[,-c(1,17)],mvnTest = "royston",univariatePlot = "qqplot")
result=mvn(new_des[,-c(1,17)],mvnTest = "royston",univariatePlot = "histogram")
result=mvn(new_des[,-c(1,17)],mvnTest = "royston",univariateTest="SW",desc = TRUE)
result$univariateNormality
result$multivariateNormality
#install.packages("Compositional")

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


# library("devtools")
# install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
# library(pairwiseAdonis)
# pairwise.adonis(iris[,1:4],iris$Species)


#**************************************************** 08/05/2018 (log de certaines variables)
attach(new_des)
plot(new_des)
names(new_des)
O2_l=log10(new_des$O2)
NO3m_l=log10(new_des$NO3_m)
NO2_l=log10(new_des$NO2)
NH4_l=log10(new_des$NH4)
#new_des$NO3         #présence de valeur nulle
#new_des$Chla        #présence de valeur nulle
#new_des$Phos       #présence de valeur nulle
#new_des$Si         # présence de val négative
################################## 09/07/2018
Fe_l=log10(new_des$Fe)
SSD_l=log10(new_des$SSD)

# des_log=data.frame(Lat,Long,T,Sal,Chla,O2_l,NO3m_l,
#                    NO3,NO2_l,NH4_l,Fe_l,SSD_l,Phos,
#                    Si,depth,month)
# result=mvn(des_log[ , -16] ,mvnTest = "royston",univariatePlot = "histogram")


################################
des_log=data.frame(Lat,Long,T,Sal,Chla,O2_l,NO3m_l,
                   NO3,NO2_l,NH4_l,Fe,SSD,Phos,
                   Si,depth)

result=mvn(des_log ,mvnTest = "royston",univariatePlot = "histogram")

names(des_log[,-16])
plot(new_des[,-1])
plot(des_log)

 
adonis(as.dist(jaccard_abundance) ~Lat+Long+T+Sal+Chla+O2_l+NO3m_l+NO3+NO2_l+NH4_l+Fe+SSD+Phos+Si+depth ,data =des_log[,-16],distance="euclidean")


mdmr.res <- mdmr(des_log,D)
summary(mdmr.res)

result=mvn(des_log[ , -16] ,mvnTest = "royston",univariatePlot = "histogram")
dataC= scale(des_log[,-16],center = TRUE,scale = TRUE) %>% as.data.frame()

# ad_scaled <- adonis(as.dist(jaccard_abundance) ~Lat+Long+T+Sal+Chla+O2_l+NO3m_l+NO3+NO2_l+NH4_l+Fe+SSD+Phos+Si+depth ,data =dataC,distance="euclidean")
# ad_scaled
# str(ad_scaled)
# ad_scaled$aov.tab[1:15,]$R2
# cbind(ad_scaled$aov.tab[1:15, ]$F.Model, ad_scaled$aov.tab[1:15, ]$R2)


#centré réduire revient à la même chose

#all.equal(ad_raw$aov.tab, ad_scaled$aov.tab)
#cbind(ad_raw$aov.tab[1:15, ]$F.Model, ad_scaled$aov.tab[1:15, ]$F.Model)
#*********************************************************

#                                   27/06/2018
#____________________________________________________________________________________________________________________
#__________________________________ Les résultats des 6 feuilles excel d'Adonis sous forme graphique ________________
library(readxl)

Fraction_Adonis <- lapply(1:6, read_sheet <- function(i) {
  read_excel("Adonis_all_/Adonis/Fraction_Adonis.xlsx",
             sheet = i, col_types = rep("text", 4))
}) %>%
  bind_rows() %>%
  filter(!is.na(Matrice)) %>%
  mutate(pval = as.numeric(`P-valeur`))

ggplot(data = Fraction_Adonis, aes(x = Matrice, y = Variables,
                                   fill = cut(pval, breaks = c(0, 0.001, 0.01, 0.05, 0.1, 1)))) +
  geom_tile() +
  facet_wrap(~Fraction, scales = "free_x") +
  scale_fill_brewer(name = "P-value", direction = -1) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
#________________________________________________________________________________________________________________________
#__________________________
#_______________________________________________________________________________________________

# ---------------- Adonis sous forme de tableau -------------------------------------
# 
# adonis_result <- adonis(as.dist(jaccard_abundance) ~Lat+Long+T+Sal+Chla+O2_l+NO3m_l+NO3+NO2_l+NH4_l+Fe+SSD+Phos+Si+depth ,data =des_log[,-16],distance="euclidean")
# 
# adonis_result$aov.tab %>% tidy() %>% mutate(Matrice = "bouh", Fraction = "5-20") #exemple


#__________________________________ Les résultats des 6 feuilles excel de bioenv sous forme graphique ________________

Fraction_Stepwise <- lapply(1:6, read_sheet <- function(i) {
  read_excel("Adonis_all_/Adonis/Fraction_Stepwise.xlsx",
             sheet = i, col_types = rep("text", 4))
}) %>%
  bind_rows() %>%
  filter(!is.na(Matrice)) %>%
  mutate(Selected = Stepwise == 1)

ggplot(data = Fraction_Stepwise, aes(x = Matrice, y = Variables,
                                     fill = Selected)) +
  geom_tile() +
  facet_wrap(~Fraction, scales = "free_x") +
  scale_fill_grey() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
#----------Create the correlation heatmap with ggplot2
#______________________________________________________________________________________
library(reshape2)
melted_cormat <- melt(jaccard_abundance)
head(melted_cormat)
library(ggplot2)
ggplot(data = melted_cormat, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile()
#__________________________________________________________________________________


names(new_des)
desi=new_des[,-c(1,17)]
names(desi)
plot(desi,col=2)
plot(des_log,col=4)

library(MVN)
result=mvn(new_des[,-c(1,17)],mvnTest = "royston",univariatePlot = "histogram")
result2=mvn(des_log,mvnTest = "royston",univariatePlot = "histogram")



#####################################################################################
################################################### 11/06/2018 #####################

#_______ RENVOIE d'un DATA FRAME des méthodes ADONIS pour une FRACTION (0-0.2)


formulas <- list(
  as.formula(as.dist(jaccard_abundance) ~Lat+Long+T+Sal+Chla+O2_l+NO3m_l+NO3+NO2_l+NH4_l+Fe+SSD+Phos+Si+depth),
  as.formula(as.dist(ochiai_abundance) ~Lat+Long+T+Sal+Chla+O2_l+NO3m_l+NO3+NO2_l+NH4_l+Fe+SSD+Phos+Si+depth),
  as.formula(as.dist(sorensen_abundance) ~Lat+Long+T+Sal+Chla+O2_l+NO3m_l+NO3+NO2_l+NH4_l+Fe+SSD+Phos+Si+depth),
  as.formula(as.dist(simka_jaccard_abundance) ~Lat+Long+T+Sal+Chla+O2_l+NO3m_l+NO3+NO2_l+NH4_l+Fe+SSD+Phos+Si+depth),
  as.formula(as.dist(chord_hellinger_prevalence) ~Lat+Long+T+Sal+Chla+O2_l+NO3m_l+NO3+NO2_l+NH4_l+Fe+SSD+Phos+Si+depth),
  as.formula(as.dist(jaccard_canberra_prevalence) ~Lat+Long+T+Sal+Chla+O2_l+NO3m_l+NO3+NO2_l+NH4_l+Fe+SSD+Phos+Si+depth),
  as.formula(as.dist(kulczynski_prevalence) ~Lat+Long+T+Sal+Chla+O2_l+NO3m_l+NO3+NO2_l+NH4_l+Fe+SSD+Phos+Si+depth),
  as.formula(as.dist(ochiai_prevalence) ~Lat+Long+T+Sal+Chla+O2_l+NO3m_l+NO3+NO2_l+NH4_l+Fe+SSD+Phos+Si+depth),
  as.formula(as.dist(whittaker_prevalence) ~Lat+Long+T+Sal+Chla+O2_l+NO3m_l+NO3+NO2_l+NH4_l+Fe+SSD+Phos+Si+depth),
  as.formula(as.dist(simka_jaccard_prevalence) ~Lat+Long+T+Sal+Chla+O2_l+NO3m_l+NO3+NO2_l+NH4_l+Fe+SSD+Phos+Si+depth),
  as.formula(as.dist(sorensen_braycurtis_prevalence) ~Lat+Long+T+Sal+Chla+O2_l+NO3m_l+NO3+NO2_l+NH4_l+Fe+SSD+Phos+Si+depth)
  
)

models <- lapply(formulas, adonis, data = des_log, distance="euclidean")

getSummary <- function(output_adonis) {
  SumsOfSqs=output_adonis$aov.tab[1:15, ]$SumsOfSqs
  F.Model=output_adonis$aov.tab[1:15, ]$F.Model
  R2=output_adonis$aov.tab[1:15, ]$R2
  P=output_adonis$aov.tab[1:15, ][[6]]
  variable <- c("Lat","Long","T","Sal","Chla","O2_l","NO3m_l","NO3",
                "NO2_l","NH4_l","Fe","SSD","Phos","Si","depth")
  tab <- as.data.frame(cbind(SumsOfSqs,F.Model,R2,P, variable))
  tab
}

summaries <- lapply(models, getSummary)
names(summaries) <- 
  c("jaccard_abundance", "ochiai_abundance", "sorensen_abundance", "simka_jaccard_abundance",
    "chord_hellinger_prevalence", "jaccard_canberra_prevalence", "kulczynski_prevalence",
    "ochiai_prevalence", "whittaker_prevalence", "simka_jaccard_prevalence","sorensen_braycurtis_prevalence")

for (k in seq_along(summaries))
  summaries[[k]]$model <- names(summaries)[[k]]

summaries_all <- Reduce("rbind", summaries)
write.csv(summaries_all, "TABLEAU.csv")










