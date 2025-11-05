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
current_directory = file.path("C:/Users/antoi/Downloads","synthetic_data_paper-main/synthetic_data_paper-main")
data_path = file.path(current_directory,"data")
setwd(data_path)
d <- 5
K <- 3
df <- read.csv("menobalance_data.csv")
df <- df[,c(1,3,4,5,7)]
varNamesDis <- NULL
varNamesCont <- c("Age", "BMI", "BP", "TotalC", "HDLC")
number_of_samples <- 1000
synthetic_data_gmcm <- generate_synthetic_data_GMCM(
  d = d,
  X = df,
  K = K,
  n_samples = number_of_samples,
  varNamesDis = varNamesDis,
  varNamesCont = varNamesCont,
  distNames = rep("bic", d),
  parametric = FALSE,
  kernel = "gaussian"
)




library(GGally)
library(ggplot2)

# Your existing code
colnames(synthetic_data_gmcm) <- colnames(df)
synthetic_data_1 <- synthetic_data_gmcm
synthetic_data_1$synthetic <- "yes"
df_1 <- df
df_1$synthetic <- "no"
df_tot <- rbind(df_1, synthetic_data_1[sample(dim(df_1)[1]),])

# Generating the ggpairs plot
p <- ggpairs(df_tot, columns = 1:5, 
             aes(color = synthetic, alphas=0.5),
             upper = list(continuous = "points"))
# Blue -> synthetic
# Red -> Real
# Display the plot
print(p)