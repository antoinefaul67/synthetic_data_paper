
library(vineclust)
library(statip)
library(mclust)
library(clue)
library(GMCM)

set.seed(51125)
# Wisconsin Breast Cancer Data
data_wisc <- read.csv("http://archive.ics.uci.edu/ml/machine-learning-databases/breast-cancer-wisconsin/wdbc.data", header = FALSE)
data_wisc_selected <- data_wisc[,c(15,27,29,30)]


d <- dim(data_wisc_selected)[2]
cdf_from_kde = function(f){
  return (function(x){integrate(f,lower = -Inf,upper = x, rel.tol=0.001)$value})
}

transform_to_uniform <- function(d, X, kernel){
  kdes = list()
  for (i in 1:d){
    kdes[[i]] = cdf_from_kde(f=densityfun(X[,i], kernel=kernel))
    uniform_X = X
    uniform_X[,i] = sapply(X[,i],kdes[[i]])
  }
  return (uniform_X)
}

U <- transform_to_uniform (d=d, X=data_wisc_selected, kernel="gaussian")
# Observed data
uhat <- Uhat(U)
# The model should be fitted multiple times using different starting estimates
start.theta <- choose.theta(uhat, m = 2) # Random starting estimate
res <- fit.full.GMCM(u = uhat, theta = start.theta,
                     method = "NM", max.ite = 3000,
                     reltol = 1e-2, trace = TRUE)

Khat <- apply(get.prob(uhat, theta = res), 1, which.max)
confusion.matrix <- table("Khat" = Khat, "K" = data_wisc$V2) # Note, some components have been swapped
ARI_GMCM <- adjustedRandIndex(Khat,data_wisc$V2)
AMI_GMCM <- AMI(Khat,data_wisc$V2)
fit_gmm <- Mclust(data_wisc_selected, G=2)
ARI_GMM <- adjustedRandIndex(fit_gmm$classification,data_wisc$V2)
AMI_GMM <-AMI(fit_gmm$classification,data_wisc$V2)
fit<- vcmm(data=data_wisc_selected, total_comp=2)
ARI_VCMM <- adjustedRandIndex(fit$cluster,data_wisc$V2)
AMI_VCMM <- AMI(fit$cluster,data_wisc$V2)
k_means_wisc <- kmeans(x= data_wisc_selected, centers=2)
ARI_Kmeans <- adjustedRandIndex(k_means_wisc$cluster,data_wisc$V2)
AMI_Kmeans <- AMI(k_means_wisc$cluster,data_wisc$V2)
confusion_matrix_k_means <- table(k_means_wisc$cluster,data_wisc$V2)
confusion.matrix.vcmm <- table(fit$cluster,data_wisc$V2)
confusion.matrix.gmm <- table(fit_gmm$classification, data_wisc$V2)

assignment_k_means <- solve_LSAP(confusion_matrix_k_means, maximum = TRUE)
confusion_rearranged_k_means <- confusion_matrix_k_means[, assignment_k_means]
assignment <- solve_LSAP(confusion.matrix, maximum = TRUE)
assignment_vcmm <- solve_LSAP(confusion.matrix.vcmm, maximum = TRUE)
assignment_gmm <- solve_LSAP(confusion.matrix.gmm, maximum = TRUE)
confusion_rearranged <- confusion.matrix[, assignment]
confusion_rearranged_vcmm <- confusion.matrix.vcmm[, assignment_vcmm]
confusion_rearranged_gmm <- confusion.matrix.gmm[, assignment_gmm]
# Misclassification rate

misclassification_rate_k_means <- 1 - sum(diag(confusion_rearranged_k_means))/sum(confusion_rearranged_k_means)
cat("Misclassification rate:", misclassification_rate_k_means, "\n")

misclassification_rate_gmcm <- 1 - sum(diag(confusion_rearranged))/sum(confusion_rearranged)
cat("Misclassification rate:", misclassification_rate, "\n")

misclassification_rate_vcmm <- 1 - sum(diag(confusion_rearranged_vcmm))/sum(confusion_rearranged_vcmm)
cat("Misclassification rate VCMM:", misclassification_rate_vcmm, "\n")

misclassification_rate_gmm <- 1 - sum(diag(confusion_rearranged_gmm))/sum(confusion_rearranged_gmm)
cat("Misclassification rate gmm:", misclassification_rate_gmm, "\n")

# Experience with simulated datasets



run_simulation <- function(n = 1000, m = 3, d = 2) {
  # Simulate data
  sim <- SimulateGMCMData(n = n, m = m, d = d)
  
  # Transform data
  X <- sim$u
  X[,1] <- qnorm(X[,1])
  X[,2] <- qgamma(X[,2], shape=2, rate=1)
  
  d=d
  U <- transform_to_uniform(d=d, X=X, kernel="gaussian")
  
  # Observed data
  uhat <- Uhat(U)
  
  # Fit GMCM model
  start.theta <- choose.theta(uhat, m = m)
  res <- fit.full.GMCM(u = uhat, theta = start.theta,
                       method = "NM", max.ite = 3000,
                       reltol = 1e-2, trace = FALSE)
  Khat <- apply(get.prob(uhat, theta = res), 1, which.max)
  confusion.matrix <- table("Khat" = Khat, "K" = sim$K)
  
  # Fit VCMM model
  fit <- vcmm(data=X, total_comp=m)
  confusion.matrix.vcmm <- table(fit$cluster, sim$K)
  
  # Fit GMM model
  fit_gmm <- Mclust(X, G=m)
  confusion.matrix.gmm <- table(fit_gmm$classification, sim$K)
  
  # Fit K-means model
  
  fit_k_means <- kmeans(x= X, centers=m)
  confusion.matrix.kmeans <- table(fit_k_means$cluster,sim$K)
  
  # Rearrange matrices
  assignment <- solve_LSAP(confusion.matrix, maximum = TRUE)
  assignment_vcmm <- solve_LSAP(confusion.matrix.vcmm, maximum = TRUE)
  assignment_gmm <- solve_LSAP(confusion.matrix.gmm, maximum = TRUE)
  assignment_k_means <- solve_LSAP(confusion.matrix.kmeans, maximum = TRUE)
  confusion_rearranged <- confusion.matrix[, assignment]
  confusion_rearranged_vcmm <- confusion.matrix.vcmm[, assignment_vcmm]
  confusion_rearranged_gmm <- confusion.matrix.gmm[, assignment_gmm]
  confusion_rearranged_k_means <- confusion.matrix.kmeans[, assignment_k_means]
  
  
  # Calculate misclassification rates
  misclassification_rate <- 1 - sum(diag(confusion_rearranged)) / sum(confusion_rearranged)
  misclassification_rate_vcmm <- 1 - sum(diag(confusion_rearranged_vcmm)) / sum(confusion_rearranged_vcmm)
  misclassification_rate_gmm <- 1 - sum(diag(confusion_rearranged_gmm)) / sum(confusion_rearranged_gmm)
  misclassification_rate_kmeans <- 1 - sum(diag(confusion_rearranged_k_means)) / sum(confusion_rearranged_k_means)
  
  #AMI
  AMI_GMCM <- AMI(Khat,sim$K)
  AMI_GMM <- AMI(fit_gmm$classification, sim$K)
  AMI_VCMM<- AMI(fit$cluster,sim$K)
  AMI_Kmeans<- AMI(fit_k_means$cluster,sim$K)
  #ARI
  ARI_GMCM <- adjustedRandIndex(Khat,sim$K)
  ARI_GMM <- adjustedRandIndex(fit_gmm$classification, sim$K)
  ARI_VCMM <- adjustedRandIndex(fit$cluster,sim$K)
  ARI_Kmeans <- adjustedRandIndex(fit_k_means$cluster,sim$K)
  print("Done")
  return(c(misclassification_rate, AMI_GMCM, ARI_GMCM, misclassification_rate_vcmm, AMI_VCMM,ARI_VCMM, misclassification_rate_gmm, AMI_GMM, ARI_GMM, misclassification_rate_kmeans, AMI_Kmeans, ARI_Kmeans))
}

n_repetitions <- 100

# Store results
misclassification_results <- replicate(n_repetitions, run_simulation())

mean_misclassification_rate_GMCM <- misclassification_results[1,]
mean_AMI_GMCM <- misclassification_results[2,]
mean_ARI_GMCM <- misclassification_results[3,]
mean_misclassification_rate_VCMM <- misclassification_results[4,]
mean_AMI_VCMM <- misclassification_results[5,]
mean_ARI_VCMM <- misclassification_results[6,]
mean_misclassification_rate_GMM <- misclassification_results[7,]
mean_AMI_GMM <- misclassification_results[8,]
mean_ARI_GMM <- misclassification_results[9,]
mean_misclassification_rate_K_means <- misclassification_results[10,]
mean_AMI_K_means <- misclassification_results[11,]
mean_ARI_K_means <- misclassification_results[12,]