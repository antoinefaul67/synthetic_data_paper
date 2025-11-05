current_directory<- file.path("C:/Users/antoi/Downloads","synthetic_data_paper-main/synthetic_data_paper-main")
data_directory <- file.path(current_directory,"data")
data_set <- "menobalance"
df = read.csv(file.path(data_directory,"menobalance_data.csv"))
set.seed(51125)
#df = df[,c(1,2,3,4,5,7)]
#df[,2] = as.factor(df[,2])
paramQuantileTrans<- function(data, varNames, distNames){
  # Initialize output list
  outList<- list()
  
  # Distributions that are currently possible
  possibleDist<- c("norm", "gamma", "lnorm")
  
  for(i in 1:length(varNames)){
    
    if(distNames[i] %in% c("aic", "bic")){
      # Choose distribution data-driven
      
      # Initialize data of performances of the distribution and list of models
      perfTab<- data.frame(dist=possibleDist, perf=NA)
      modelDist<- list()
      
      
      for(j in possibleDist){
        if(j=="norm" | (j  %in% c("gamma", "lnorm") & min(data[[varNames[i]]])>0)){
          # Fit distribution model
          modelDist[[j]]<- fitdistrplus::fitdist(data[, varNames[i]], distr=j, 
                                                 method="mle")
          
          # Extract performance (AIC or BIC)
          perfTab$perf[perfTab$dist==j]<- modelDist[[j]][[distNames[i]]]
        }
      }
      perfTab<- na.omit(perfTab)
      
      # Select distribution with best performance
      chosenDist<- perfTab$dist[perfTab$perf<=min(perfTab$perf)]
      chosenDist<- chosenDist[1]
      
      estimatesVar<- modelDist[[chosenDist]]$estimate
      
    } else{
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

## transform data to standard normal distribution
quantileTransData<- function(model, data){
  varNames<- names(model)
  for(i in varNames){
    distName<- model[[i]]$dist
    distEst<- model[[i]]$estimates
    if(distName=="norm"){
      data[, i]<- qnorm(pnorm(data[, i], distEst[1], distEst[2]))
    } else if(distName=="gamma"){
      data[, i]<- qnorm(pgamma(data[, i], distEst[1], rate = distEst[2]))
    } else if(distName=="lnorm"){
      data[, i]<- qnorm(plnorm(data[, i], distEst[1], distEst[2]))
    }
  }
  return(data)
}

## re-transform data to original scale
quantileRetransData<- function(model, data){
  varNames<- names(model)
  for(i in varNames){
    distName<- model[[i]]$dist
    distEst<- model[[i]]$estimates
    if(i %in% colnames(data)){
      if(distName=="norm"){
        data[, i]<- qnorm(pnorm(data[, i]), distEst[1], distEst[2])
      } else if(distName=="gamma"){
        data[, i]<- qgamma(pnorm(data[, i]), distEst[1], rate = distEst[2])
      } else if(distName=="lnorm"){
        data[, i]<- qlnorm(pnorm(data[, i]), distEst[1], distEst[2])
      }
    }
  }
  return(data)
}


condMultiNorm<- function(respNames, pred, varNamesTrans, distNamesTrans, data){
  # Define number of responses and number of variables (responses and predictors)
  nResp<- length(respNames)
  nVar<- nResp+length(pred)
  varNamesTrans <- varNamesTrans
  distNamesTrans<- distNamesTrans
  # Define training data
  dataTrain<- data[, c(respNames, pred)]
  modelTrans<- paramQuantileTrans(dataTrain, varNamesTrans, distNamesTrans)
  dataTrain<- quantileTransData(modelTrans, dataTrain)
  # Compute mean vector and covariance matrix
  mean_vector = apply(dataTrain, 2, mean)
  cov_matrix = cov(dataTrain)
  
  # Compute covariance matrix of response conditional on predictor values
  # (independent of actual predictor values!)
  if (length(respNames)==1){
    cov_cond = cov_matrix[1:nResp, 1:nResp] - as.matrix(t(cov_matrix[1:nResp, (nResp+1):nVar])) %*%
      as.matrix(solve(cov_matrix[(nResp+1):nVar, (nResp+1):nVar])) %*% 
      as.matrix(cov_matrix[(nResp+1):nVar, 1:nResp])
  }
  else{
    cov_cond = cov_matrix[1:nResp, 1:nResp] - as.matrix(cov_matrix[1:nResp, (nResp+1):nVar]) %*%
      as.matrix(solve(cov_matrix[(nResp+1):nVar, (nResp+1):nVar])) %*% 
      as.matrix(cov_matrix[(nResp+1):nVar, 1:nResp])
  }
  outList<- list(input=list(type="condNv",respNames=respNames, pred=pred), 
                 model=list(means=mean_vector, covAll=cov_matrix, cov=cov_cond, modelTrans=modelTrans))
}
## prediction
predMeanCondMultiNorm<- function(model, data){
  # Extract info from model
  modelTrans <-  model$model$modelTrans
  data <- quantileTransData(modelTrans, data)
  respNames<- model$input$respNames
  predNames<- model$input$pred
  mean_vector<- model$model$means
  cov_matrix<- model$model$covAll
  
  # Define number of responses and number of variables (responses and predictors)
  nResp<- length(respNames)
  nVar<- nResp+length(predNames)
  
  # Define data used for prediction
  
  
  dataPred<- data[, c(respNames, predNames)]
  
  # predicted mean
  for(i in 1:nrow(dataPred)){
    dataSub<- dataPred[i, ]
    if (length(respNames)==1){
      meanCondSub<- mean_vector[1:nResp] + as.matrix(t(cov_matrix[1:nResp, (nResp+1):nVar])) %*%
        as.matrix(solve(cov_matrix[(nResp+1):nVar, (nResp+1):nVar])) %*%
        t(as.matrix(dataSub[(nResp+1):nVar] - mean_vector[(nResp+1):nVar]))
    }
    else{
      meanCondSub<- mean_vector[1:nResp] + as.matrix(cov_matrix[1:nResp, (nResp+1):nVar]) %*%
        as.matrix(solve(cov_matrix[(nResp+1):nVar, (nResp+1):nVar])) %*%
        t(as.matrix(dataSub[(nResp+1):nVar] - mean_vector[(nResp+1):nVar]))
    }
    
    meanCondSub<- t(meanCondSub)
    
    if(i==1){
      meanCond<- meanCondSub
    } else{
      meanCond<- rbind(meanCond, meanCondSub)
    }
  }
  
  # return predicted (e.g., conditional) mean
  return(meanCond)
}

# Simulate data from predicted multivariate normal distribution
predSimDist<- function(model, dataSub, nSim=NULL, seed=131820){
  # load library
  library(MASS)
  modelTrans <- model$model$modelTrans
  # Initialize output list
  outList<- list()
  dataSim <- list()
  # Extract info from model
  respNames<- model$input$respNames
  predNames<- model$input$pred
  
  # Gaussian Copula
  if(model$input$type=="condNv"){
    # Extract covariance matrix
    covMat<- model$model$cov
    dataSub<- dataSub[, c(respNames, predNames)]
    predMeans<- predMeanCondMultiNorm(model, dataSub)
  }
  
  # Simulate data based on predictions
  for (i in 1:(dim(predMeans)[1])){
    dataSim[[i]]<- mvrnorm(nSim, predMeans[i,], covMat)
    dataSim[[i]]<- quantileRetransData(model=modelTrans, data=dataSim[[i]])
  }
  outList$mean<- predMeans
  outList$sim<- dataSim
  
  return(outList)
}

score2Risk<- function(vec){
  age = as.numeric(vec[1])
  smoker = as.numeric(vec[2])
  systBP = as.numeric(vec[3])
  totalChol = as.numeric(vec[4])
  hdlChol = as.numeric(vec[5])
  diabetes = as.numeric(vec[6])
  
  cage = (age-60)/5
  csbp = (systBP-120)/20
  ctchol = (totalChol-6)/1
  chdl = (hdlChol-1.3)/0.5
  
  # female scale for low risk region
  scale1 = -0.7380
  scale2 = 0.7019
  
  beta = as.matrix(c(0.4648, 0.7744, 0.3131, 0.8096, 0.1002, -0.2606, -0.1088, 
                     -0.0277, -0.0226, 0.0613, -0.1272))
  
  var = as.matrix(c(cage, smoker, csbp, diabetes, ctchol, chdl, cage*smoker, 
                    cage*csbp, cage*ctchol, cage*chdl, cage*diabetes))
  
  L = t(beta) %*% var
  Prob = 1 - 0.9776^exp(L)
  
  #calibrating probability
  Prob_calibrated = 1-exp(-exp(scale1 + scale2*log(-log(1-Prob))))
  return(Prob_calibrated*100)
}
library(scoringRules)

compute_performance <- function(dataTest, sim){
  nSim <- dim(sim[[1]])[1]
  crpsPred <- rep(0,length(sim))
  logsPred<- rep(0,length(sim))
  wcrpsPred <- rep(0,length(sim))
  for (i in 1:length(sim)){
    sim[[i]] = as.data.frame(sim[[i]])
    sim[[i]]$Age<- rep(dataTest[i,]$Age, nSim)
    sim[[i]]$Smoker<- rep(dataTest[i,]$Smoker, nSim)
    sim[[i]]$Diabetes<- rep(dataTest[i,]$Diabetes, nSim)
    sim[[i]]<-  sim[[i]][, c("Age", "Smoker", "BP", "TotalC", "HDLC", "Diabetes")]
    sim[[i]]$Score2<- apply(sim[[i]],1,score2Risk)
    est <- sim[[i]]$Score2
    true <- dataTest$Score2[i]
    crpsPred[i] <- scoringRules::crps_sample(y=true,dat=est)
    logsPred[i]<- scoringRules::logs_sample(y=true, dat=est)
    wcrpsPred[i] <- scoringRules::twcrps_sample(y=true, dat=est, a=5)
  }
  return (list(crpsPred=crpsPred,logsPred=logsPred, wcrpsPred = wcrpsPred))
}
library(dplyr)

# Define a function to compute the mean CRPS and log score for multiple splits and datasets
compute_mean_metrics <- function(df, synthetic_folder_paths, train_ratio = 0.8, n_splits = 10, varNamesTrans) {
  result_list <- list()
  
  # Iterate over synthetic data folders
  for (folder_path in synthetic_folder_paths) {
    synthetic_files <- list.files(folder_path, pattern = "\\.csv$|\\.rds$", full.names = TRUE)
    
    folder_metrics <- sapply(synthetic_files, function(synthetic_file) {
      if (grepl("\\.csv$", synthetic_file)) {
        synthetic_data <- read.csv(synthetic_file)
      } else if (grepl("\\.rds$", synthetic_file)) {
        synthetic_data <- readRDS(synthetic_file)
      }
      
      metrics_values <- replicate(n_splits, {
        # Create train-test split
        train_index <- sample(seq_len(nrow(df)), floor(train_ratio * nrow(df)), replace = FALSE)
        
        dataTrainReel <- df[train_index, ]
        dataTest <- df[-train_index, ]
        
        dataTrainSynthetic <- synthetic_data
        
        # Fit models
        model_GC_synth <- condMultiNorm(respNames = c("BP", "TotalC", "HDLC"),
                                        pred = c("Age", "BMI"),
                                        varNamesTrans =varNamesTrans,
                                        distNamesTrans = rep("aic",5),
                                        data = dataTrainSynthetic)
        
        # Make predictions
        pred_synth <- predSimDist(model_GC_synth, dataSub = dataTest, nSim = 1000)
        dataTest$Score2 <- apply(dataTest, 1, score2Risk)
        
        # Compute performance
        perf_synth <- compute_performance(dataTest = dataTest, sim = pred_synth$sim)
        
        # Calculate the mean CRPS and log score for synthetic data
        c(mean_crps = mean(perf_synth$crpsPred), mean_logs = mean(perf_synth$logsPred), mean_wrcps = (perf_synth$wcrpsPred) )
      })
      
      # Return means of metrics across all splits
      rowMeans(metrics_values)
    })
    print(paste0("done",folder_path))
    # Store the result for the folder
    result_list[[basename(folder_path)]] <- folder_metrics
  }
  
  return(result_list)
}

# Usage
varNamesTrans = c("BP","TotalC","HDLC","Age","BMI")
distNamesTrans = rep("aic",5)
GMCM_path <- file.path(current_directory, data_set,"GMCM")
marginal_path <- file.path(current_directory, data_set,"marginals")
GC_path <-  file.path(current_directory, data_set,"GC synthetizer")
CTGAN_path <- file.path(current_directory, data_set,"CTGAN synthetizer")
Vines_path <- file.path(current_directory, data_set,"Vines")
synthetic_folder_paths <- c(GC_path, CTGAN_path,GMCM_path, marginal_path, Vines_path)
#df<-df[,c(1,2,3,4,5,7)]
mean_results <-compute_mean_metrics(df, synthetic_folder_paths, train_ratio = 0.8, n_splits = 10,varNamesTrans=varNamesTrans)

mean(mean_results$`GC synthetizer`["mean_crps",])

# Performance from the real data
num_splits <- 100 # You can adjust this as needed

# Initialize a vector to store CRPS values for each split
crps_values <- numeric(num_splits)
logs_values <- numeric(num_splits)
wcrps_values <- numeric(num_splits)
# Loop over the number of splits
for (i in 1:num_splits) {
  
  # Split the data into training and test sets
  set.seed(i) # Set seed for reproducibility
  train_indices <- sample(seq_len(nrow(df)), size = 0.8 * nrow(df))
  dataTrainReel <- df[train_indices, ]
  dataTest <- df[-train_indices, ]
  
  # Train the model on the training data
  model_GC_real <- condMultiNorm(respNames = c("BP", "TotalC", "HDLC"),
                                 pred = c("Age", "BMI"),varNamesTrans=c("BP", "TotalC", "HDLC","Age","BMI"),distNamesTrans=rep("aic",5),
                                 data = dataTrainReel)
  
  # Generate predictions on the test data
  pred_real <- predSimDist(model_GC_real, dataSub = dataTest, nSim = 1000)
  dataTest$Score2 <- apply(dataTest,1,score2Risk)
  # Compute performance
  perf_real <- compute_performance(dataTest = dataTest, sim = pred_real$sim)
  
  # Store the CRPS for this split
  crps_values[i] <- mean(perf_real$crpsPred)
  logs_values[i] <- mean(perf_real$logsPred)
  wcrps_values[i] <- mean(perf_real$wcrpsPred)
}

# Compute the overall mean CRPS across all splits
mean_crps_real <- mean(crps_values)
mean_logs_real <- mean(logs_values)
mean_wrcps_real <- mean(wcrps_values)