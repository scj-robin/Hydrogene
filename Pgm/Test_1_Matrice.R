library(vegan)
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(gdata)
library(readxl)
library(mice)
library(RCurl)

set.seed(1234)
data.wd <- "/home/diabakhate/HydroGen_Abdou"
setwd(dir=data.wd) 
#--
data1= read_excel("~/HydroGen_Abdou/Genocenoses_env_parameters_all_tara_initial.xlsx")
#data1= read_excel("~/HydroGen_Abdou/Genocenoses_env_parameters_all_tara.xls")

names(data1)
#unique(data1$Genocenose)

data.5_20 = data1[which((as.character(data1$Fraction) == "43952") == TRUE),]     #dim(77,20)
size_fraction4 = "43952" ########c'est le 5-20
rownames(data.5_20) = data.5_20$Station

design = data.5_20

import_data(size_fraction4, samples = rownames(design))
design <- design[metagenomic_sample, ]  #nouveau jeu de taille (77,20)
dim(design)
summary(design)
#unique(design$Genocenose)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Lecture des covariables
Genocenoses_env_parameters_all_tara <- data1
CovarTot = as.data.frame(Genocenoses_env_parameters_all_tara);
CovarTot$Temp = CovarTot$T; CovarTot$T = c()
#------------------------
# Transformation log
CovarTot$Sal = log10(CovarTot$Sal); CovarTot$Chla = log10(1+CovarTot$Chla); CovarTot$O2 = log10(CovarTot$O2);
CovarTot$NO3_m = log10(CovarTot$NO3_m); CovarTot$NO3 = log10(1+CovarTot$NO3); CovarTot$NO2 = log10(CovarTot$NO2);
CovarTot$NH4 = log10(CovarTot$NH4); CovarTot$Fe = log10(CovarTot$Fe); CovarTot$Phos = log10(1+CovarTot$Phos);
CovarTot$Si = log10(CovarTot$Si);

# #------------------------
# # Imputation
CovarImput = mice::complete(mice(CovarTot, method="norm.predict", m=5))
fraction = '43952'           #5-20
Covar = CovarImput[which(CovarImput$Fraction==fraction), ] # on s'intéresse à cette sous data frame de fraction
names(Covar)
#Covar = Covar[, c(3:13,14:26,27,30)] #Notre jeu de donnée avec les colonnes qui nous intéresse
Covar = Covar[, c(3,4,5,6:16,17,20)]
n = nrow(Covar); p = ncol(Covar)
head(Covar)
#summary(Covar)

Group = c(4,1,1,rep(3,11),2,2) # 1 = geographie, 2 = physique, 3= chimie,4=Genocenose
#Group = c(rep(4,11),1,1,rep(3,11),2,2)
GroupName = c('geographie','physique', 'chimie','Genocenose')
G = max(Group)    
GroupOrder = order(Group)         # on ordonne le groupe suivant les chiffres 1-2-3
GroupOrdered = Group[order(Group)]
CovarOrdered = Covar[order(Group)] ##Ordonne le data frame (les variables) suivant 'geographie', 'physique', 'chimie'
head(CovarOrdered)

 CovarOrdered$Genocenose=as.factor(CovarOrdered$Genocenose)
# str(CovarOrdered)
 unique(CovarOrdered$Genocenose)  #les modalités concernées

########################  On filtre, variance non nulle
# Var = apply(CovarOrdered, 2, var)
# CovarOrderedReduced = CovarOrdered[, which(Var > 0)]
# head(CovarOrderedReduced)

# Group = c(1,1,2,2,rep(3,11),rep(4,6))
# GroupName = c('geographie','physique', 'chimie','Genocenose')


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
library(parallel)
nb_cores <- 7

#
formulas <- list(as.formula(as.dist(jaccard_abundance) ~.))

models <- mclapply(formulas, adonis, data = CovarOrdered, mc.cores = nb_cores)
# models <- mclapply(formulas, adonis, data = CovarOrderedReduced, mc.cores = nb_cores)

getSummary <- function(output_adonis) {
  
  TotalSumsOfSqs = output_adonis$aov.tab$SumsOfSqs[length(output_adonis$aov.tab$SumsOfSqs)]
  ResidSumsOfSqs = output_adonis$aov.tab$SumsOfSqs[length(output_adonis$aov.tab$SumsOfSqs)-1]
  ResidMeanSqs = output_adonis$aov.tab$MeanSqs[length(output_adonis$aov.tab$MeanSqs)-1]
  ResidDf = output_adonis$aov.tab$Df[length(output_adonis$aov.tab$Df)-1]
  Group = data.frame()
  GroupMeanSqs = GroupSumsOfSqs = GroupF = GroupDf = GroupR2 = rep(0, G)
  for (g in 1:G){
    GroupList = which(GroupOrdered==g)
    # GroupDf[g] = length(GroupList)
    
    GroupDf[g] = sum(output_adonis$aov.tab$Df[GroupList])
    GroupSumsOfSqs[g] = sum(output_adonis$aov.tab$SumsOfSqs[GroupList])
    GroupMeanSqs[g] = GroupSumsOfSqs[g] / GroupDf[g]
    GroupF[g] = GroupMeanSqs[g] / ResidMeanSqs                  #look 
    # GroupF[g] = (ResidDf*output_adonis$aov.tab$SumsOfSqs[g])/(ResidSumsOfSqs*output_adonis$aov.tab$Df[g])
    GroupR2[g] = GroupSumsOfSqs[g] / TotalSumsOfSqs
  }
  
  # Permutations
  B = 999; GroupPval = rep(0, G)   
  for (g in 1:G){
    Fpermut = rep(0, B)
    for (b in 1:B){
      if(b%%round(sqrt(B))==0){cat(b, '')}
      GroupList = which(GroupOrdered==g)
      Permut = order(runif(n))
      CovarPermut = CovarOrdered
      CovarPermut[, GroupList] = CovarPermut[Permut, which(GroupOrdered==g)]
      ADONISPermut = adonis(formula(output_adonis), CovarPermut, permutations=2)
      ResidSumsOfSqsPermut = ADONISPermut$aov.tab$SumsOfSqs[length(ADONISPermut$aov.tab$SumsOfSqs)-1]
      GroupSumsOfSqsPermut = sum(ADONISPermut$aov.tab$SumsOfSqs[GroupList])
      Fpermut[b] = GroupSumsOfSqsPermut / length(GroupList) / (ResidSumsOfSqsPermut / ResidDf)
    }
     hist(Fpermut, breaks=sqrt(B), main=paste(GroupName[g], ': F=', round(Fpermut[g], 3))); abline(v=Fpermut[g], col=2);
    GroupPval[g] = (1+sum(Fpermut>=Fpermut[g]))/(1+B)
  }
  
  GroupAov = as.data.frame(cbind(GroupDf, GroupSumsOfSqs, GroupMeanSqs, GroupF, GroupR2, GroupPval, GroupName))
  row.names(GroupAov) = GroupName
  GroupAov
}
summaries <- mclapply(models, getSummary, mc.cores = nb_cores)
summaries

