import_data <- function(size_fraction, samples = NULL){
oldwd <- getwd()  
  
if(size_fraction=="0.8-5"){
  setwd(dir="/home/diabakhate/HydroGen_Abdou/metagenomic_covariates_Le Bon/Donnees/simka_matrices_17-05-22/Simka_0.8-5") 
}

if(size_fraction=="43952"){  # 5-20
  setwd(dir="/home/diabakhate/HydroGen_Abdou/metagenomic_covariates_Le Bon/Donnees/simka_matrices_17-05-22/Simka_5-20") 
}

if(size_fraction=="180-2000"){
  setwd(dir="/home/diabakhate/HydroGen_Abdou/metagenomic_covariates_Le Bon/Donnees/simka_matrices_17-05-22/Simka_180-2000") 
}

if(size_fraction=="0.22-3"){
  setwd(dir="/home/diabakhate/HydroGen_Abdou/metagenomic_covariates_Le Bon/Donnees/simka_matrices_17-05-22/Simka_0.22-3") 
}

if(size_fraction=="0-0.2"){
  setwd(dir="/home/diabakhate/HydroGen_Abdou/metagenomic_covariates_Le Bon/Donnees/simka_matrices_17-05-22/Simka_0-0.2") 
}

if(size_fraction=="20-180"){
  setwd(dir="/home/diabakhate/HydroGen_Abdou/metagenomic_covariates_Le Bon/Donnees/simka_matrices_17-05-22/Simka_20-180") 
}

############################################################ Import data 

if(size_fraction=="0.8-5" | size_fraction=="43952" | size_fraction=="180-2000"){
  jaccard_abundance2 = read.table(file="mat_abundance_jaccard.csv",sep="")
  if(dim(jaccard_abundance2)[2]==1){
    jaccard_abundance2 = read.table(file="mat_abundance_jaccard.csv",sep=";")
    sep = ";"
  }else{
    sep = ""
  }
  delete=c(which(duplicated(jaccard_abundance2)==TRUE)) # remove duplicates
  jaccard_abundance2 = jaccard_abundance2[-delete,-delete]
  metagenomic_sample <<- as.character(jaccard_abundance2[-1,1])
  
  jaccard_abundance2 = read.table(file="mat_abundance_jaccard.csv",sep=sep)
  ab_jaccard_abundance2 = read.table(file="mat_abundance_ab-jaccard.csv",sep=sep)
  braycurtis_abundance2 = read.table(file="mat_abundance_braycurtis.csv",sep=sep)
  ab_ochiai_abundance2 = read.table(file="mat_abundance_ab-ochiai.csv",sep=sep)
  ab_sorensen_abundance2 = read.table(file="mat_abundance_ab-sorensen.csv",sep=sep)
  simka_jaccard_abundance2 = read.table(file="mat_abundance_simka-jaccard.csv",sep=sep)
  
  chord_prevalence2 = read.table(file="mat_presenceAbsence_chord.csv",sep=sep)
  jaccard_prevalence2 = read.table(file="mat_presenceAbsence_jaccard.csv",sep=sep)
  kulczynski_prevalence2 = read.table(file="mat_presenceAbsence_kulczynski.csv",sep=sep)
  ochiai_prevalence2 = read.table(file="mat_presenceAbsence_ochiai.csv",sep=sep)
  whittaker_prevalence2 = read.table(file="mat_presenceAbsence_whittaker.csv",sep=sep)
  simka_jaccard_prevalence2 = read.table(file="mat_presenceAbsence_simka-jaccard.csv",sep=sep)
  sorensen_braycurtis_prevalence2 = read.table(file="mat_presenceAbsence_sorensen-braycurtis.csv",sep="\t") 
  
  ############################################################ Distance matrices
  
  jaccard_abundance2 = unname(as.matrix(jaccard_abundance2)[c(-1,-delete),c(-1,-delete)])
  ab_jaccard_abundance2 = unname(as.matrix(ab_jaccard_abundance2)[c(-1,-delete),c(-1,-delete)])
  braycurtis_abundance2 = unname(as.matrix(braycurtis_abundance2)[c(-1,-delete),c(-1,-delete)])
  ab_ochiai_abundance2 = unname(as.matrix(ab_ochiai_abundance2)[c(-1,-delete),c(-1,-delete)])
  ab_sorensen_abundance2 = unname(as.matrix(ab_sorensen_abundance2)[c(-1,-delete),c(-1,-delete)])
  simka_jaccard_abundance2 = unname(as.matrix(simka_jaccard_abundance2)[c(-1,-delete),c(-1,-delete)])
  
  
  chord_prevalence2 = unname(as.matrix(chord_prevalence2)[c(-1,-delete),c(-1,-delete)])
  jaccard_prevalence2 = unname(as.matrix(jaccard_prevalence2)[c(-1,-delete),c(-1,-delete)])
  kulczynski_prevalence2 = unname(as.matrix(kulczynski_prevalence2)[c(-1,-delete),c(-1,-delete)])
  ochiai_prevalence2 = unname(as.matrix(ochiai_prevalence2)[c(-1,-delete),c(-1,-delete)])
  whittaker_prevalence2 = unname(as.matrix(whittaker_prevalence2)[c(-1,-delete),c(-1,-delete)])
  simka_jaccard_prevalence2 = unname(as.matrix(simka_jaccard_prevalence2)[c(-1,-delete),c(-1,-delete)])
  sorensen_braycurtis_prevalence2 = unname(as.matrix(sorensen_braycurtis_prevalence2)[c(-1,-delete),c(-1,-delete)])
  
  n = dim(jaccard_abundance2)[1]
  
  jaccard_abundance <<- matrix(NA,nrow=n,ncol=n)
  ab_jaccard_abundance <<- matrix(NA,nrow=n,ncol=n)
  braycurtis_abundance <<- matrix(NA,nrow=n,ncol=n)
  ab_ochiai_abundance <<- matrix(NA,nrow=n,ncol=n)
  ab_sorensen_abundance <<- matrix(NA,nrow=n,ncol=n)
  simka_jaccard_abundance <<- matrix(NA,nrow=n,ncol=n)
  
  chord_prevalence <<- matrix(NA,nrow=n,ncol=n)
  jaccard_prevalence <<- matrix(NA,nrow=n,ncol=n)
  kulczynski_prevalence <<- matrix(NA,nrow=n,ncol=n)
  ochiai_prevalence <<- matrix(NA,nrow=n,ncol=n)
  whittaker_prevalence <<- matrix(NA,nrow=n,ncol=n)
  simka_jaccard_prevalence <<- matrix(NA,nrow=n,ncol=n)
  sorensen_braycurtis_prevalence <<- matrix(NA,nrow=n,ncol=n)
  
  for(j in 1:n){
    jaccard_abundance[,j] <<- as.numeric(jaccard_abundance2[,j])
    ab_jaccard_abundance[,j] <<- as.numeric(ab_jaccard_abundance2[,j])
    braycurtis_abundance[,j] <<- as.numeric(braycurtis_abundance2[,j])
    ab_ochiai_abundance[,j] <<- as.numeric(ab_ochiai_abundance2[,j])
    ab_sorensen_abundance[,j] <<- as.numeric(ab_sorensen_abundance2[,j])
    simka_jaccard_abundance[,j] <<- as.numeric(simka_jaccard_abundance2[,j])
    
    chord_prevalence[,j] <<- as.numeric(chord_prevalence2[,j])
    jaccard_prevalence[,j] <<- as.numeric(jaccard_prevalence2[,j])
    kulczynski_prevalence[,j] <<- as.numeric(kulczynski_prevalence2[,j])
    ochiai_prevalence[,j] <<- as.numeric(ochiai_prevalence2[,j])
    whittaker_prevalence[,j] <<- as.numeric(whittaker_prevalence2[,j])
    simka_jaccard_prevalence[,j] <<- as.numeric(simka_jaccard_prevalence2[,j])
    sorensen_braycurtis_prevalence[,j] <<- as.numeric(sorensen_braycurtis_prevalence2[,j])
  }
  
  ############################################################ Name rows and columns correctly
  
  matrices.list <<- c("jaccard_abundance", "ab_jaccard_abundance", "braycurtis_abundance", 
                     "ab_ochiai_abundance", "ab_sorensen_abundance", "simka_jaccard_abundance", 
                     "chord_prevalence", "jaccard_prevalence", "kulczynski_prevalence", 
                     "ochiai_prevalence", "whittaker_prevalence", "simka_jaccard_prevalence",
                     "sorensen_braycurtis_prevalence")
  for (mat in matrices.list) {
    dist.mat <- get(mat)
    rownames(dist.mat) <- colnames(dist.mat) <- metagenomic_sample
    if (!is.null(samples)) {
      complete_sample <- intersect(metagenomic_sample, samples)
      dist.mat <- dist.mat[complete_sample, complete_sample]
    }
    assign(x = mat, value = dist.mat, envir = .GlobalEnv)
  }
  
  if (!is.null(samples)) {
    metagenomic_sample <<- complete_sample
    if (length(metagenomic_sample) == 0) {
      warning("No sample selected, check the value of `samples` argument")
    }
  }
  
}else{
  sep = ""
  jaccard_abundance2 = read.table(file="mat_abundance_jaccard.csv",sep=sep)
  ochiai_abundance2 = read.table(file="mat_abundance_ochiai.csv",sep=sep)
  sorensen_abundance2 = read.table(file="mat_abundance_sorensen.csv",sep=sep)
  simka_jaccard_abundance2 = read.table(file="mat_abundance_simka-jaccard.csv",sep=sep)
  
  chord_hellinger_prevalence2 = read.table(file="mat_presenceAbsence_chord-hellinger.csv",sep=sep)
  jaccard_canberra_prevalence2 = read.table(file="mat_presenceAbsence_jaccard-canberra.csv",sep=sep)
  kulczynski_prevalence2 = read.table(file="mat_presenceAbsence_kulczynski.csv",sep=sep)
  ochiai_prevalence2 = read.table(file="mat_presenceAbsence_ochiai.csv",sep=sep)
  whittaker_prevalence2 = read.table(file="mat_presenceAbsence_whittaker.csv",sep=sep)
  simka_jaccard_prevalence2 = read.table(file="mat_presenceAbsence_simka-jaccard.csv",sep=sep)
  sorensen_braycurtis_prevalence2 = read.table(file="mat_presenceAbsence_sorensen-braycurtis.csv",sep=sep)
  
  metagenomic_sample <<- row.names(jaccard_abundance2) #### ICI ON s'intéresse à la première matrice
  
  ############################################################ Distance matrices
  
  jaccard_abundance2 = unname(as.matrix(jaccard_abundance2))
  ochiai_abundance2 = unname(as.matrix(ochiai_abundance2))
  sorensen_abundance2 = unname(as.matrix(sorensen_abundance2))
  simka_jaccard_abundance2 = unname(as.matrix(simka_jaccard_abundance2))
  
  
  chord_hellinger_prevalence2 = unname(as.matrix(chord_hellinger_prevalence2))
  jaccard_canberra_prevalence2 = unname(as.matrix(jaccard_canberra_prevalence2))
  kulczynski_prevalence2 = unname(as.matrix(kulczynski_prevalence2))
  ochiai_prevalence2 = unname(as.matrix(ochiai_prevalence2))
  whittaker_prevalence2 = unname(as.matrix(whittaker_prevalence2))
  simka_jaccard_prevalence2 = unname(as.matrix(simka_jaccard_prevalence2))
  sorensen_braycurtis_prevalence2 = unname(as.matrix(sorensen_braycurtis_prevalence2))
  
  n = dim(jaccard_abundance2)[1]
  
  jaccard_abundance <<- matrix(NA,nrow=n,ncol=n)
  ochiai_abundance <<- matrix(NA,nrow=n,ncol=n)
  sorensen_abundance <<- matrix(NA,nrow=n,ncol=n)
  simka_jaccard_abundance <<- matrix(NA,nrow=n,ncol=n)
  
  
  chord_hellinger_prevalence <<- matrix(NA,nrow=n,ncol=n)
  jaccard_canberra_prevalence <<- matrix(NA,nrow=n,ncol=n)
  kulczynski_prevalence <<- matrix(NA,nrow=n,ncol=n)
  ochiai_prevalence <<- matrix(NA,nrow=n,ncol=n)
  whittaker_prevalence <<- matrix(NA,nrow=n,ncol=n)
  simka_jaccard_prevalence <<- matrix(NA,nrow=n,ncol=n)
  sorensen_braycurtis_prevalence <<- matrix(NA,nrow=n,ncol=n)
  
  for(j in 1:n){
    jaccard_abundance[,j] <<- as.numeric(jaccard_abundance2[,j])
    ochiai_abundance[,j] <<- as.numeric(ochiai_abundance2[,j])
    sorensen_abundance[,j] <<- as.numeric(sorensen_abundance2[,j])
    simka_jaccard_abundance[,j] <<- as.numeric(simka_jaccard_abundance2[,j])
    
    chord_hellinger_prevalence[,j] <<- as.numeric(chord_hellinger_prevalence2[,j])
    jaccard_canberra_prevalence[,j] <<- as.numeric(jaccard_canberra_prevalence2[,j])
    kulczynski_prevalence[,j] <<- as.numeric(kulczynski_prevalence2[,j])
    ochiai_prevalence[,j] <<- as.numeric(ochiai_prevalence2[,j])
    whittaker_prevalence[,j] <<- as.numeric(whittaker_prevalence2[,j])
    simka_jaccard_prevalence[,j] <<- as.numeric(simka_jaccard_prevalence2[,j])
    sorensen_braycurtis_prevalence[,j] <<- as.numeric(sorensen_braycurtis_prevalence2[,j])
  }
  
  ############################################################ Name rows and columns correctly
  
  matrices.list <<- c("jaccard_abundance", "ochiai_abundance", "sorensen_abundance", 
                      "simka_jaccard_abundance", "chord_hellinger_prevalence", 
                     "jaccard_canberra_prevalence", "kulczynski_prevalence", 
                     "ochiai_prevalence", "whittaker_prevalence", "simka_jaccard_prevalence",
                     "sorensen_braycurtis_prevalence")
  
  for (mat in matrices.list) {
    dist.mat <- get(mat)
    rownames(dist.mat) <- colnames(dist.mat) <- metagenomic_sample
    if (!is.null(samples)) {
      complete_sample <- intersect(metagenomic_sample, samples)
      dist.mat <- dist.mat[complete_sample, complete_sample]
    }
    assign(x = mat, value = dist.mat, envir = .GlobalEnv)
  }
  
  if (!is.null(samples)) {
    metagenomic_sample <<- complete_sample
    if (length(metagenomic_sample) == 0) {
      warning("No sample selected, check the value of `samples` argument")
    }
  }
}

############################################################ Restore working directory

setwd(oldwd)

}

spectral.clustering.new <- function(A, normalised = TRUE, score = FALSE, K = 2, adj = FALSE){
  
  ### preparation 
  n = dim(A)[1]
  iso.A = isolate(A)
  iso.seq = iso.A$isolate
  noniso.seq = iso.A$nonisolate
  A.noniso = A[noniso.seq, noniso.seq]
  labels = rep(0, n)
  
  ### svd
  n = dim(A.noniso)[1]
  eig.A = eigen(A.noniso)
  if(score == F){
    U = matrix(0, nrow = n, ncol = K)	
  }
  if(score == T){
    U = matrix(0, nrow = n, ncol = (K-1))
  }
  
  ### get U matrix
  if(adj == F & score == F){
    L = laplacian(A = A.noniso, normalised = normalised)
    eig = eigen(L, symmetric = TRUE)
    
    for(i in 1:K){
      U[, i] = eig$vectors[, (n + 1 - i)]
    }
  }
  
  if(adj == T & score == F){
    ordered.vec = eig.A$vectors[, order(abs(eig.A$values), decreasing = T)]
    for(i in 1:K){
      U[, i] = ordered.vec[, (n + 1 - i)]
    }
  }
  
  if(score == T){
    ordered.vec = eig.A$vectors[, order(abs(eig.A$values), decreasing = T)]
    benchmark = ordered.vec[,1] + 1e-5
    for(i in 2:K){
      U[, (i-1)] = eig.A$vectors[, order(abs(eig.A$values), decreasing = F)][, i]/benchmark
    }
  }
  
  ### k-means
  U = scale(U, center = F)
  temp = unique(U, margin = 2)
  if(dim(temp)[1] < K){stop('FAIL!')}
  k.means = kmeans(U, centers = K, nstart = 1000, iter.max = 20)
  labels[noniso.seq] = k.means$cluster
  
  return(labels)
}


