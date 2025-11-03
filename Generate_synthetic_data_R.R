# Organized generation of synthetic data

# Necessary functions.

library(fitdistrplus)
library(GGally)
library(GMCM)
library(statip)

set.seed(031125)
cdf_from_kde = function(f){
  return (function(x){integrate(f,lower = -Inf,upper = x, rel.tol=0.001)$value})
}
inverse_function<- function (q,f, lower, upper) {
  uniroot(function(x){f(x) - q}, lower=lower, upper=upper)$root
}

paramQuantileTrans<- function(data, varNames, distNames){
  # Initialize output list
  outList<- list()
  possibleDist <- c("norm", "gamma", "lnorm","cauchy","weibull")
  
  for(i in 1:length(varNames)) {
    
    if(distNames[i] %in% c("aic", "bic")) {
      # Choose distribution data-driven
      
      # Initialize data of performances of the distribution and list of models
      perfTab <- data.frame(dist = possibleDist, perf = NA)
      modelDist <- list()
      
      for(j in possibleDist) {
        if((j %in% c("norm","cauchy","weibull"))| (j %in% c("gamma", "lnorm") & min(data[[varNames[i]]]) > 0)) {
          
          # Fit distribution model with error handling
          fit_result <- tryCatch({
            fitdistrplus::fitdist(data[, varNames[i]], distr = j, method = "mle")
          }, error = function(e) {
            return(NULL)
          })
          
          if(!is.null(fit_result)) {
            # Extract performance (AIC or BIC)
            modelDist[[j]] <- fit_result
            perfTab$perf[perfTab$dist == j] <- fit_result[[distNames[i]]]
          }
        }
      }
      
      perfTab <- na.omit(perfTab)
      
      # Select distribution with best performance if available
      if(nrow(perfTab) > 0) {
        chosenDist <- perfTab$dist[perfTab$perf <= min(perfTab$perf)]
        chosenDist <- chosenDist[1]
        
        estimatesVar <- modelDist[[chosenDist]]$estimate
      } else {
        # Handle case where no valid distribution is found (optional)
        warning("No valid distribution found for variable: ", varNames[i])
      }
    }
    else{
      # Fit model
      chosenDist<- distNames[i]
      estimatesVar<- fitdistrplus::fitdist(data[, varNames[i]], distr=distNames[i], 
                                           method="mle")$estimate
    }
    
    # save it to outputlist
    outList[[varNames[i]]]<- list(dist=chosenDist, estimates=estimatesVar)
  }
  
  return(outList)
}

UnifTransData <- function(model, data){
  varNames<- names(model)
  for(i in varNames){
    distName<- model[[i]]$dist
    distEst<- model[[i]]$estimates
    if(i %in% colnames(data)){
      if(distName=="norm"){
        data[, i]<- pnorm(data[, i], distEst[1], distEst[2])
      } else if(distName=="gamma"){
        data[, i]<- pgamma(data[, i], distEst[1], rate = distEst[2])
      } else if(distName=="lnorm"){
        data[, i]<- plnorm(data[, i], distEst[1], distEst[2])
      } else if(distName=="weibull"){
        data[, i]<- qnorm(pweibull(data[, i], distEst[1], scale = distEst[2]))
      } else if(distName=="cauchy"){
        data[, i]<- qnorm(pcauchy(data[, i], distEst[1], distEst[2]))
      }
      
    }
  }
  return(data)
}

quantileRetransData<- function(model, data){
  varNames<- names(model)
  for(i in varNames){
    distName<- model[[i]]$dist
    distEst<- model[[i]]$estimates
    if(i %in% colnames(data)){
      if(distName=="norm"){
        data[, i]<- qnorm(data[, i], distEst[1], distEst[2])
      } else if(distName=="gamma"){
        data[, i]<- qgamma(data[, i], distEst[1], rate = distEst[2])
      } else if(distName=="lnorm"){
        data[, i]<- qlnorm(data[, i], distEst[1], distEst[2])
      } else if(distName=="weibull"){
        data[, i]<- qweibull(pnorm(data[, i]), distEst[1], shape=distEst[2])
      } else if(distName=="cauchy"){
        data[, i]<- qcauchy(pnorm(data[, i]), distEst[1], distEst[2])
      }
    }
  }
  return(data)
}

generate_synthetic_data_GMCM <- function(d, X, K, n_samples, varNamesCont, distNames, varNamesDis, parametric, kernel){
  n <- nrow(X)
  kdes <- list()
  if (parametric){
    modelTrans <- paramQuantileTrans(data=X, varNames=varNamesCont, distNames=distNames)
    unif_data <- UnifTransData(model = modelTrans, data = X)
  }
  else{
    unif_data <- X
    for (i in varNamesCont){
      kdes[[i]] = cdf_from_kde(f=densityfun(X[,i], kernel=kernel))
      unif_data[,i] = sapply(X[,i],kdes[[i]])
    }
  }
  intervals <- list()
  for (v in varNamesDis) {
    # Ensure factor
    X[[v]] <- as.factor(X[[v]])
    
    # Get the counts (table) and proportions
    tbl <- table(X[[v]])
    props <- prop.table(tbl)
    
    # Cumulative proportions define the intervals' break points
    cumsum_props <- c(0, cumsum(as.numeric(props)))
    
    levels_v <- names(props)
    
    # Build list of intervals for each level
    intervals[[v]] <- data.frame(
      level = levels_v,
      left = head(cumsum_props, -1),
      right = tail(cumsum_props, -1)
    )
    interval <- intervals[[v]]
    # Create empty vector to store the sampled values
    sampled_uniform <- numeric(nrow(X))
    
    # For each observation, sample
    for (i in seq_len(nrow(X))) {
      fac <- X[[v]][i]
      inter <- interval[interval$level == fac, ]
      sampled_uniform[i] <- runif(1, min = inter$left, max = inter$right)
    }
    
    # Optionally, store in new column
    unif_data[[v]] <- sampled_uniform
  }
  #fit<- capture.output({fit.full.GMCM(Uhat(unif_data), m=K, max.ite = 1000, method = "NM", theta = choose.theta(Uhat(unif_data), m = K))})
  fit <- fit.full.GMCM(Uhat(unif_data), m=K, max.ite = 1000, method = "NM", theta = choose.theta(Uhat(unif_data), m = K))
  ## Generate uniform data following the fitted copula
  synthetic_data_uniform<-SimulateGMCMData(n=n_samples, theta = fit)$u
  synthetic_data_uniform = as.data.frame(synthetic_data_uniform)
  colnames(synthetic_data_uniform) <- colnames(unif_data)
  reTransData <- synthetic_data_uniform
  if (parametric){
    reTransData <- quantileRetransData(model= modelTrans, data=synthetic_data_uniform)
  }
  else{
    for (i in varNamesCont){
      reTransData[,i] <- sapply(synthetic_data_uniform[,i],inverse_function, f= kdes[[i]], lower = min(X[,i])-0.5*mean(X[,i]), upper=max(X[,i])+0.5*mean(X[,i]))
    }
  }
  for (v in varNamesDis){
    interval <- intervals[[v]]
    breaks <- c(interval$left[1], interval$right)
    levels <- interval$level
    
    # Assign intervals
    reTransData[[v]] <-cut(reTransData[[v]],
                           breaks = breaks,
                           labels = levels,
                           include.lowest = TRUE,
                           right = TRUE)
    reTransData[[v]] <- as.numeric(as.character(reTransData[[v]]))
  }
  return (reTransData)
}
generate_synthetic_data_vines<- function(d, X, n_samples, varNamesCont, distNames, varNamesDis, parametric, kernel){
  n <- nrow(X)
  kdes <- list()
  if (parametric){
    modelTrans <- paramQuantileTrans(data=X, varNames=varNamesCont, distNames=distNames)
    unif_data <- UnifTransData(model = modelTrans, data = X)
  }
  else{
    unif_data <- X
    for (i in varNamesCont){
      kdes[[i]] = cdf_from_kde(f=densityfun(X[,i], kernel=kernel))
      unif_data[,i] = sapply(X[,i],kdes[[i]])
    }
  }
  intervals <- list()
  for (v in varNamesDis) {
    # Ensure factor
    X[[v]] <- as.factor(X[[v]])
    
    # Get the counts (table) and proportions
    tbl <- table(X[[v]])
    props <- prop.table(tbl)
    
    # Cumulative proportions define the intervals' break points
    cumsum_props <- c(0, cumsum(as.numeric(props)))
    
    levels_v <- names(props)
    
    # Build list of intervals for each level
    intervals[[v]] <- data.frame(
      level = levels_v,
      left = head(cumsum_props, -1),
      right = tail(cumsum_props, -1)
    )
    interval <- intervals[[v]]
    # Create empty vector to store the sampled values
    sampled_uniform <- numeric(nrow(X))
    
    # For each observation, sample
    for (i in seq_len(nrow(X))) {
      fac <- X[[v]][i]
      inter <- interval[interval$level == fac, ]
      sampled_uniform[i] <- runif(1, min = inter$left, max = inter$right)
    }
    
    # Optionally, store in new column
    unif_data[[v]] <- sampled_uniform
  }
  fit <- VineCopula::RVineStructureSelect(Uhat(unif_data), familyset = c(1,2))  # Fit R-vine with families 1 (Gaussian), 2 (t)
  #simulated <- RVineSim(1000, fit)
  #fit<- capture.output({fit.full.GMCM(Uhat(unif_data), m=K, max.ite = 1000, method = "NM", theta = choose.theta(Uhat(unif_data), m = K))})
  #fit <- vinecop(data=unif_data)
  ## Generate uniform data following the fitted copula
  #synthetic_data_uniform<-SimulateGMCMData(n=n_samples, theta = fit)$u
  synthetic_data_uniform<-VineCopula::RVineSim(n_samples, fit)
  synthetic_data_uniform = as.data.frame(synthetic_data_uniform)
  colnames(synthetic_data_uniform) <- colnames(unif_data)
  reTransData <- synthetic_data_uniform
  if (parametric){
    reTransData <- quantileRetransData(model= modelTrans, data=synthetic_data_uniform)
  }
  else{
    for (i in varNamesCont){
      reTransData[,i] <- sapply(synthetic_data_uniform[,i],inverse_function, f= kdes[[i]], lower = min(X[,i])-0.5*mean(X[,i]), upper=max(X[,i])+0.5*mean(X[,i]))
    }
  }
  for (v in varNamesDis){
    interval <- intervals[[v]]
    breaks <- c(interval$left[1], interval$right)
    levels <- interval$level
    
    # Assign intervals
    reTransData[[v]] <-cut(reTransData[[v]],
                           breaks = breaks,
                           labels = levels,
                           include.lowest = TRUE,
                           right = TRUE)
    reTransData[[v]] <- as.numeric(as.character(reTransData[[v]]))
  }
  return (reTransData)
}

generate_synthetic_data_marginals <- function(df,n_samples,kernel,contVar, disVar, distNames = rep("aic",length(contVar)+length(disVar)),kde=TRUE){
  kdes = list()
  synthetic_data <- data.frame(matrix(ncol = 0, nrow = n_samples))
  if (kde){
    for (v in contVar){
      kdes[[v]] = cdf_from_kde(f=densityfun(df[[v]], kernel=kernel))
      values <- rep(0,n_samples)
      for (i in 1:n_samples){
        values[i] <- inverse_function(q=runif(1),f= kdes[[v]], lower = min(df[[v]])-0.5*mean(df[[v]]), upper=max(df[[v]])+0.5*mean(df[[v]]))
      }
      synthetic_data[[v]]<-values
    }
  }
  else{
    for (v in contVar){
      outList<- list()
      possibleDist <- c("norm", "gamma", "lnorm","cauchy","weibull")
      perfTab <- data.frame(dist = possibleDist, perf = NA)
      modelDist <- list()
      for(j in possibleDist) {
        if((j %in% c("norm","cauchy","weibull"))| (j %in% c("gamma", "lnorm") & min(df[, v]) > 0)) {
          
          # Fit distribution model with error handling
          fit_result_1 <- tryCatch({
            fitdistrplus::fitdist(df[, v], distr = j, method = "mle")
          }, error = function(e) {
            return(NULL)
          })
          
          if(!is.null(fit_result_1)) {
            # Extract performance (AIC or BIC)
            modelDist[[j]] <- fit_result_1
            perfTab$perf[perfTab$dist == j] <- fit_result_1[["aic"]]
          }
        }
      }
      
      perfTab <- na.omit(perfTab)
      
      # Select distribution with best performance if available
      if(nrow(perfTab) > 0) {
        chosenDist <- perfTab$dist[perfTab$perf <= min(perfTab$perf)]
        chosenDist <- chosenDist[1]
        estimatesVar <- modelDist[[chosenDist]]$estimate
        
      }
      possibleDist<- c("norm", "gamma", "lnorm")
      distEst<- modelDist[[chosenDist]]$estimate
      values <- rep(0,n_samples)
      for (i in 1:n_samples){
        if(chosenDist=="norm"){
          values[i]<- qnorm(runif(1), distEst[1], distEst[2])
        } else if(chosenDist=="gamma"){
          values[i]<- qgamma(runif(1), distEst[1], rate = distEst[2])
        } else if(chosenDist=="lnorm"){
          values[i]<- qlnorm(runif(1), distEst[1], distEst[2])
        } else if(chosenDist=="weibull"){
          values[i]<- qweibull(runif(1), distEst[1], shape=distEst[2])
        } else if(chosenDist=="cauchy"){
          values[i]<- qcauchy(runif(1), distEst[1], distEst[2])
        }
      }
      synthetic_data[[v]]<-values
    }
  } 
  
  for (v in disVar){
    my_fac <- as.factor(df[[v]])
    props <- prop.table(table(my_fac))
    x = data.frame(props)
    samples <- sample(
      x$my_fac, 
      size = n_samples, 
      replace = TRUE, 
      prob = x$Freq
    )
    
    # Result: sampled vector with proportions similar to original
    synthetic_data[[v]] <- samples
  }
  synthetic_data <- synthetic_data[,colnames(df)]
  return (synthetic_data)
}


current_directory = getwd()

data_path = file.path(current_directory,"data")

setwd(data_path)

if (dataset == "menobalance") {
  df <- read.csv("menobalance_data.csv")
  df_selected <- df[, c(1,2,3,4,5,6,7)]
  varNamesDis <- c("Smoker", "Diabetes")
  varNamesCont <- c("Age", "BMI", "BP", "TotalC", "HDLC")
  d <- 7
  K <- 3
  X <- df_selected
} else if (dataset == "cleveland") {
  data_clev <- read.csv("https://archive.ics.uci.edu/ml/machine-learning-databases/heart-disease/processed.cleveland.data", header = FALSE)
  names(data_clev) <- c("age", "sex", "cp", "trestbps", "chol",
                        "fbs", "restecg", "thalach", "exang", "oldpeak",
                        "slope", "ca", "thal", "num")
  data_clev_selected <- subset(data_clev, select = -c(num))
  data_clev_selected <- subset(data_clev_selected, select = c("trestbps", "thalach", "chol", "age"))
  d <- 4
  K <- 2
  varNamesCont <- c("trestbps", "thalach", "chol", "age")
  varNamesDis <- NULL
  X <- data_clev_selected
} else if (dataset == "wisconsin") {
  data_wisc <- read.csv("http://archive.ics.uci.edu/ml/machine-learning-databases/breast-cancer-wisconsin/wdbc.data", header = FALSE)
  data_wisc_selected <- data_wisc[, c(15, 27, 29, 30)]
  d <- 4
  varNamesCont <- c("V15", "V27", "V29", "V30")
  varNamesDis <- NULL
  K <- 2
  X <- data_wisc_selected
} else if (dataset == "simulated_gmcm") {
  df <- read.csv("simulated_GMCM_data.csv")
  df_selected <- df[, 2:3]
  d <- 2
  K <- 2
  varNamesDis <- NULL
  varNamesCont <- c("X1", "X2")
  X <- df_selected
} else {
  stop("Unknown dataset")
}


number_of_synthetic_datasets <- 100
number_of_samples <- 1000
successful_runs <- 0
attempt <- 1

experiment_directory = file.path(current_directory, dataset)
subdirectories = c("Vines", "marginals", "GMCM")

# Create the subdirectories if they don't already exist
for (subdir in subdirectories) {
  dir_path = file.path(experiment_directory, subdir)
  if (!dir.exists(dir_path)) {
    dir.create(dir_path, recursive = TRUE)
  }
}

while (successful_runs < number_of_synthetic_datasets) {
  try({
    synthetic_data_vines<-  generate_synthetic_data_vines(
      d = d,
      X = X,
      n_samples = number_of_samples,
      varNamesDis = varNamesDis,
      varNamesCont = varNamesCont,
      distNames = rep("bic", d),
      parametric = FALSE,
      kernel = "gaussian"
    )
    synthetic_data_gmcm <- generate_synthetic_data_GMCM(
      d = d,
      X = X,
      K = K,
      n_samples = number_of_samples,
      varNamesDis = varNamesDis,
      varNamesCont = varNamesCont,
      distNames = rep("bic", d),
      parametric = FALSE,
      kernel = "gaussian"
    )
    synthetic_data_marginal <- generate_synthetic_data_marginals(df=X,n_samples=number_of_samples
                                                                    ,kernel="gaussian", contVar= varNamesCont,
                                                                    disVar=varNamesDis,kde=TRUE)
    colnames(synthetic_data_gmcm) <- colnames(X)
    colnames(synthetic_data_marginal) <- colnames(X)
    colnames(synthetic_data_vines) <- colnames(X)
    saveRDS(object = synthetic_data_marginal,
            file =file.path(experiment_directory,"marginals",paste0("synthetic_data_", successful_runs + 1, ".rds")))
    saveRDS(object = synthetic_data_gmcm
            , file = file.path(experiment_directory,"GMCM",paste0("synthetic_data_", successful_runs + 1, ".rds")))
    saveRDS(object = synthetic_data_vines
            , file = file.path(experiment_directory,"Vines",paste0("synthetic_data_", successful_runs + 1, ".rds")))
   
    folder_name <- "Vines"
    file_path <- file.path(folder_name, file_name)
    saveRDS(object = synthetic_data_vines, file = file_path)
    successful_runs <- successful_runs + 1
  }, silent = TRUE)
  
  attempt <- attempt + 1
  print(successful_runs)
}