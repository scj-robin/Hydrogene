rm(list=ls()); par(mfrow=c(1, 1))
library(vegan); library(dplyr); library(tidyr); library(ggplot2); library(stringr); library(gdata); library(readxl); 
library(mice)

#------------------------
# Parms
data.dir <- "../Data/"
distance.name = 'mat_abundance_jaccard'
fraction = '0.8-5'
analysis.name = paste0('frac', fraction, '_', distance.name)

#------------------------
# Lecture des distances
mat_abundance_jaccard <- read.table(paste0(data.dir, distance.name, ".csv"), sep='\t', header=T)
row.names(mat_abundance_jaccard) = mat_abundance_jaccard[, 1]; mat_abundance_jaccard = mat_abundance_jaccard[, -1]
Dist = as.dist(mat_abundance_jaccard)
dim(as.matrix(Dist))

#------------------------
# Lecture des covariables
# Genocenoses_env_parameters_all_tara <- read_excel(paste0(data.dir, "Genocenoses_env_parameters_all_tara.xls"))
# CovarTot = as.data.frame(Genocenoses_env_parameters_all_tara);
# CovarTot$Temp = CovarTot$T; CovarTot$T = c()
# #------------------------
# # Transformation log
# CovarTot$Sal = log10(CovarTot$Sal); CovarTot$Chla = log10(1+CovarTot$Chla); CovarTot$O2 = log10(CovarTot$O2);
# CovarTot$NO3_m = log10(CovarTot$NO3_m); CovarTot$NO3 = log10(1+CovarTot$NO3); CovarTot$NO2 = log10(CovarTot$NO2);
# CovarTot$NH4 = log10(CovarTot$NH4); CovarTot$Fe = log10(CovarTot$Fe); CovarTot$Phos = log10(1+CovarTot$Phos);
# CovarTot$Si = log10(CovarTot$Si);
# summary(CovarTot)
# #------------------------
# # Imputation
# CovarImput = complete(mice(CovarTot, method="norm.predict", m=5))
# summary(CovarImput)
# save(CovarImput, file=paste0(data.dir, 'Genocenoses_Imput.Rdata'))
load(paste0(data.dir, 'Genocenoses_Imput.Rdata'))

#------------------------
# Selection des stations
Covar = CovarImput[which(CovarImput$Fraction==fraction), ]
Covar = Covar[, c(4, 5, 20, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17)]
n = nrow(Covar); p = ncol(Covar)
head(Covar)

#------------------------
# Adonis
# ADONIS = adonis(Dist ~ ., Covar)

#------------------------
# Définition des groupes des variables
names(Covar)
Group = c(1, 1, 2, rep(3, 11), 2) # 1 = geographie, 2 = physique, 3 = chimie
GroupName = c('geograph', 'physique', 'chimie')
G = max(Group)
GroupOrder = order(Group)
GroupOrdered = Group[order(Group)]
CovarOrdered = Covar[order(Group)]

# Adonis
ADONIS = adonis(Dist ~ ., CovarOrdered)

# Statistiques des groupes
TotalSumsOfSqs = ADONIS$aov.tab$SumsOfSqs[length(ADONIS$aov.tab$SumsOfSqs)]
ResidSumsOfSqs = ADONIS$aov.tab$SumsOfSqs[length(ADONIS$aov.tab$SumsOfSqs)-1]
ResidMeanSqs = ADONIS$aov.tab$MeanSqs[length(ADONIS$aov.tab$MeanSqs)-1]
ResidDf = ADONIS$aov.tab$Df[length(ADONIS$aov.tab$Df)-1]
Group = data.frame()
GroupMeanSqs = GroupSumsOfSqs = GroupF = GroupDf = GroupR2 = rep(0, G)
for (g in 1:G){
   GroupList = which(GroupOrdered==g)
   GroupDf[g] = length(GroupList)
   GroupSumsOfSqs[g] = sum(ADONIS$aov.tab$SumsOfSqs[GroupList])
   GroupMeanSqs[g] = GroupSumsOfSqs[g] / GroupDf[g]
   GroupF[g] = GroupMeanSqs[g] / ResidMeanSqs
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
      ADONISPermut = adonis(Dist ~ ., CovarPermut, permutations=2)
      ResidSumsOfSqsPermut = ADONISPermut$aov.tab$SumsOfSqs[length(ADONISPermut$aov.tab$SumsOfSqs)-1]
      GroupSumsOfSqsPermut = sum(ADONISPermut$aov.tab$SumsOfSqs[GroupList])
      Fpermut[b] = GroupSumsOfSqsPermut / length(GroupList) / (ResidSumsOfSqsPermut / ResidDf)
      }
   hist(Fpermut, breaks=sqrt(B), main=paste(GroupName[g], ': F=', round(Fpermut[g], 3))); abline(v=Fpermut[g], col=2); 
   GroupPval[g] = (1+sum(Fpermut>=Fpermut[g]))/(1+B)
}
GroupAov = as.data.frame(cbind(GroupDf, GroupSumsOfSqs, GroupMeanSqs, GroupF, GroupR2, GroupPval))
row.names(GroupAov) = GroupName
GroupAov


############################
# TODO
# Creer une fonction qui prenne en entree Dist, Covar, Group (B=999, )