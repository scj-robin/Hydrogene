library(vegan)
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(gdata)
library(readxl)

data.wd <- "C:/Users/pc_user/Documents/HydroGen_Abdou"
setwd(dir=data.wd) 
#--
data1= read_excel("~/HydroGen_Abdou/Genocenoses_env_parameters_all_tara.xlsx")
names(data1)


data.0_0.2 = data1[which((as.character(data1$Fraction) == "0-0.2") == TRUE),]    #dim(70,20)
data.0.22_3 = data1[which((as.character(data1$Fraction) == "0.22-3") == TRUE),]  #dim(105,20)
data.0.8_5 = data1[which((as.character(data1$Fraction) == "0.8-5") == TRUE),]    #dim(131,20)
data.5_20 = data1[which((as.character(data1$Fraction) == "43952") == TRUE),]     #dim(77,20)
data.20_180 = data1[which((as.character(data1$Fraction) == "20-180") == TRUE),]  #dim(126,20)
data.180_2000 = data1[which((as.character(data1$Fraction) == "180-2000") == TRUE),] #dim(135,20)

rownames(data.0_0.2) = data.0_0.2$Station
rownames(data.0.22_3) = data.0.22_3$Station
rownames(data.0.8_5) = data.0.8_5$Station
rownames(data.5_20) = data.5_20$Station
rownames(data.20_180) = data.20_180$Station
rownames(data.180_2000) = data.180_2000$Station

#**************************************
size_fraction1 = "0-0.2"
size_fraction2 = "0.22-3"
size_fraction3 = "0.8-5"
size_fraction4 = "43952" ########c'est le 5-20
size_fraction5 = "20-180"
size_fraction6 = "180-2000"

design = data.5_20

import_data(size_fraction4, samples = rownames(design))

design <- design[metagenomic_sample, ]  #nouveau jeu de taille (70,20)
dim(design)
#------------------------

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#------------------------
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
summary(CovarTot)
# #------------------------
# # Imputation
library(mice)
CovarImput = complete(mice(CovarTot, method="norm.predict", m=5))
summary(CovarImput)
#save(CovarImput, file=paste0(data.dir, 'Genocenoses_Imput.Rdata'))
#load(paste0(data.dir, 'Genocenoses_Imput.Rdata'))

#------------------------
# Selection des stations
fraction = '43952'
Covar = CovarImput[which(CovarImput$Fraction==fraction), ] # on s'intéresse à cette sous data frame de fraction
names(Covar)
Covar = Covar[, c(3,6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16,17,20,4, 5)] #Notre jeu de donnée avec les colonnes qui nous intéresse
n = nrow(Covar); p = ncol(Covar)
head(Covar)


#------------------------
# Définition des groupes des variables
#Group = c(3,4, 4, 1, rep(2, 11), 1) # 4 = geographie, 1 = physique, 2= chimie,3=Genocenose
Group = c(1,rep(2,11))
#GroupName = c('physique', 'chimie','geographie','Genocenose')
GroupName = c('Genocenose','Chimie')
G = max(Group)    # la valeur 3
GroupOrder = order(Group)         # on ordonne le groupe suivant les chiffres 1-2-3
GroupOrdered = Group[order(Group)]
CovarOrdered = Covar[order(Group)] ##Ordonne le data frame (les variables) suivant 'geographie', 'physique', 'chimie'
head(CovarOrdered)

# Adonis
# Dist = as.dist(jaccard_abundance)              # Pour une matrice donnée
# ADONIS = adonis(Dist ~ ., CovarOrdered)
# ADONIS
#
# # Statistiques des groupes
# TotalSumsOfSqs = ADONIS$aov.tab$SumsOfSqs[length(ADONIS$aov.tab$SumsOfSqs)]
# ResidSumsOfSqs = ADONIS$aov.tab$SumsOfSqs[length(ADONIS$aov.tab$SumsOfSqs)-1]
# ResidMeanSqs = ADONIS$aov.tab$MeanSqs[length(ADONIS$aov.tab$MeanSqs)-1]
# ResidDf = ADONIS$aov.tab$Df[length(ADONIS$aov.tab$Df)-1]
# Group = data.frame()
# GroupMeanSqs = GroupSumsOfSqs = GroupF = GroupDf = GroupR2 = rep(0, G)
# for (g in 1:G){
#   GroupList = which(GroupOrdered==g)
#   GroupDf[g] = length(GroupList)
#   GroupSumsOfSqs[g] = sum(ADONIS$aov.tab$SumsOfSqs[GroupList])
#   GroupMeanSqs[g] = GroupSumsOfSqs[g] / GroupDf[g]
#   GroupF[g] = GroupMeanSqs[g] / ResidMeanSqs
#   GroupR2[g] = GroupSumsOfSqs[g] / TotalSumsOfSqs
# }
#
# # Permutations
# B = 999; GroupPval = rep(0, G)  #B=999
# for (g in 1:G){
#   Fpermut = rep(0, B)
#   for (b in 1:B){
#     if(b%%round(sqrt(B))==0){cat(b, '')}
#     GroupList = which(GroupOrdered==g)
#     Permut = order(runif(n))
#     CovarPermut = CovarOrdered
#     CovarPermut[, GroupList] = CovarPermut[Permut, which(GroupOrdered==g)]
#     ADONISPermut = adonis(Dist ~ ., CovarPermut, permutations=2)
#     ResidSumsOfSqsPermut = ADONISPermut$aov.tab$SumsOfSqs[length(ADONISPermut$aov.tab$SumsOfSqs)-1]
#     GroupSumsOfSqsPermut = sum(ADONISPermut$aov.tab$SumsOfSqs[GroupList])
#     Fpermut[b] = GroupSumsOfSqsPermut / length(GroupList) / (ResidSumsOfSqsPermut / ResidDf)
#   }
#   hist(Fpermut, breaks=sqrt(B), main=paste(GroupName[g], ': F=', round(Fpermut[g], 3))); abline(v=Fpermut[g], col=2);
#   GroupPval[g] = (1+sum(Fpermut>=Fpermut[g]))/(1+B)
# }
# GroupAov = as.data.frame(cbind(GroupDf, GroupSumsOfSqs, GroupMeanSqs, GroupF, GroupR2, GroupPval))
# row.names(GroupAov) = GroupName
# GroupAov


#_________________________________________________________________________________________________________
#_________________________________________________________________________________________________________
########################################## 18/07/2018
############################
# Les résultats de adonis pour les variables concaténées
library(parallel)
#nb_cores <- 7

#
formulas <- list(
as.formula(as.dist(jaccard_abundance) ~.),
as.formula(as.dist(ab_jaccard_abundance) ~.),
as.formula(as.dist(braycurtis_abundance) ~.),
as.formula(as.dist(ab_ochiai_abundance) ~.),
as.formula(as.dist(ab_sorensen_abundance) ~.),
as.formula(as.dist(simka_jaccard_abundance) ~.),
as.formula(as.dist(chord_prevalence) ~.),
as.formula(as.dist(jaccard_prevalence) ~.),
as.formula(as.dist(kulczynski_prevalence) ~.),
as.formula(as.dist(ochiai_prevalence) ~.),
as.formula(as.dist(whittaker_prevalence) ~.),
as.formula(as.dist(simka_jaccard_prevalence) ~.),
as.formula(as.dist(sorensen_braycurtis_prevalence) ~.))

# formulas <- list(
#   as.formula(as.dist(jaccard_abundance) ~.),
#   as.formula(as.dist(ochiai_abundance) ~.),
#   as.formula(as.dist(sorensen_abundance) ~.),
#   as.formula(as.dist(simka_jaccard_abundance) ~.),
#   as.formula(as.dist(chord_hellinger_prevalence) ~.),
#   as.formula(as.dist(jaccard_canberra_prevalence) ~.),
#   as.formula(as.dist(kulczynski_prevalence) ~.),
#   as.formula(as.dist(ochiai_prevalence) ~.),
#   as.formula(as.dist(whittaker_prevalence) ~.),
#   as.formula(as.dist(simka_jaccard_prevalence) ~.),
#   as.formula(as.dist(sorensen_braycurtis_prevalence) ~.))

models <- mclapply(formulas, adonis, data = CovarOrdered)
#models <- mclapply(formulas, adonis, data = CovarOrdered, mc.cores = nb_cores)
getSummary <- function(output_adonis) {

  TotalSumsOfSqs = output_adonis$aov.tab$SumsOfSqs[length(output_adonis$aov.tab$SumsOfSqs)]
  ResidSumsOfSqs = output_adonis$aov.tab$SumsOfSqs[length(output_adonis$aov.tab$SumsOfSqs)-1]
  ResidMeanSqs = output_adonis$aov.tab$MeanSqs[length(output_adonis$aov.tab$MeanSqs)-1]
  ResidDf = output_adonis$aov.tab$Df[length(output_adonis$aov.tab$Df)-1]
  Group = data.frame()
  GroupMeanSqs = GroupSumsOfSqs = GroupF = GroupDf = GroupR2 = rep(0, G)
  for (g in 1:G){
    GroupList = which(GroupOrdered==g)
    GroupDf[g] = length(GroupList)
    GroupSumsOfSqs[g] = sum(output_adonis$aov.tab$SumsOfSqs[GroupList])
    GroupMeanSqs[g] = GroupSumsOfSqs[g] / GroupDf[g]
    GroupF[g] = GroupMeanSqs[g] / ResidMeanSqs
    GroupR2[g] = GroupSumsOfSqs[g] / TotalSumsOfSqs
  }

  # Permutations
  B = 999; GroupPval = rep(0, G)   # NORMALEMENT B=999 MAIS ÇA PREND DU TEMPS À COMPILER C POURQUOI J'AI CHANGÉ LA VALEUR
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

summaries <- mclapply(models, getSummary)
#summaries <- mclapply(models, getSummary, mc.cores = nb_cores)

#
names(summaries) <-
  c("jaccard_abundance","ab_jaccard_abundance","braycurtis_abundance","ab_ochiai_abundance","ab_sorensen_abundance",
    "simka_jaccard_abundance","chord_prevalence","jaccard_prevalence","kulczynski_prevalence","ochiai_prevalence",
    "whittaker_prevalence","simka_jaccard_prevalence","sorensen_braycurtis_prevalence")

# names(summaries) <- c("jaccard_abundance","ochiai_abundance","sorensen_abundance","simka_jaccard_abundance","chord_hellinger_prevalence",
#                       "jaccard_canberra_prevalence","kulczynski_prevalence","ochiai_prevalence","whittaker_prevalence",
#                       "simka_jaccard_prevalence","sorensen_braycurtis_prevalence")

for (k in seq_along(summaries))
  summaries[[k]]$model <- names(summaries)[[k]]

summaries_all <- Reduce("rbind", summaries)
rownames(summaries_all) <- 1:nrow(summaries_all)

head(summaries_all)

write.table(summaries_all,"tab5-20.csv",sep = ";")

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#############################################################################################

#________________________ Les résultats des 6 feuilles excel d'Adonis sous forme graphique ____
 library(readxl)
#
Adonis_Pull <- lapply(1:6, read_sheet <- function(i) {
 # read_excel("Adonis_all_/Adonis/Adonis_Pull.xlsx",
  #           sheet = i, col_types = rep("text", 4))
  read_excel("C:/Users/pc_user/Documents/HydroGen_Abdou/AdonisStat/Genocenose_Chimie.xlsx",
                       sheet = i, col_types = rep("text", 4))
}) %>%
  bind_rows() %>%
  filter(!is.na(Matrice)) %>%
  mutate(pval = as.numeric(`GroupPval`))

ggplot(data = Adonis_Pull, aes(x = Matrice, y = GroupName,
                                   fill = cut(pval, breaks = c(0, 0.001, 0.01, 0.05, 0.1, 1)))) +
  geom_tile() +
  facet_wrap(~Fraction, scales = "free_x") +
  scale_fill_brewer(name = "P-value", direction = -1) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
