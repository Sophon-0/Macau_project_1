
###################################################################### Sparse Group Lasso #########################################################################
###################################################################################################################################################################
library(SGL)
# For type="linear" should be a list with $x$ an input matrix of dimension n-obs
# by p-vars, and $y$ a length $n$ response vector.
# For type="cox" should be a list with x as before, time, an n-vector
# corresponding to failure/censor times, and status, an n-vector indicating failure (1) or censoring (0).

train_SGL <- function(group_index, train_feat, train_obs, xTrain_feat, xTrain_obs, alpha) {
  
  train_feat <- as.matrix(train_feat)
  train_obs <- as.matrix(train_obs)
  xTrain_feat <- as.matrix(xTrain_feat)
  data = list(x = train_feat, y = train_obs)
  
  # train SGL model 
  fit <- SGL(data, group_index, type = "linear", alpha=alpha, nlam = 20)
  
  xTrain_pred <- predictSGL(fit, xTrain_feat ) ## rows are samples, columns are each step of determination of beta coefficients 
   
  xTrain_cIDX <- c()
  for (i in 1:ncol(fit$beta)) {
    xTrain_cIDX <- c(xTrain_cIDX, cIDX(xTrain_pred[,i], xTrain_obs))  ## Harrell c-index or concordance C
    ## fit$beta is a matrix of the beta weights across the different steps for determining the betas
    ## we save all the concordance index for every step of beta for cross validation set
  }
  
  lambda_num <- as.numeric(which(xTrain_cIDX == max(xTrain_cIDX)))
  ## choose the lambda for the highest c index in the cross validation set
  
  return(list(model=fit,
              lambda_num=lambda_num,
              xTrain_cIDX=xTrain_cIDX))
}


RUN_SGL <- function(index,group_index,fold, iteration, alpha) {
  res <- RES[index, ]
  res <- res[!is.na(res)]
  if(length(res) < 30) { fold<-3 }
  
  pearson_table <- c()
  RMSE_table <- c()
  for (i in 1:iteration) {    
    crossVal <- crossvalidation(names(res), nFold=fold, seed=F)
    observation <- c()
    prediction <- c()
    
    for (i in 1:fold) {
      
      train_obs <- res[as.character(crossVal[[i]]$train)]
      train_feat <- as.matrix(features[names(train_obs), ])
      rownames(train_feat) <- names(train_obs)
      train_feat[is.na(train_feat)] <- 0  ## replace NA by 0
      
      # prepare cross-trainset
      xTrain_obs <- res[as.character(crossVal[[i]]$xTrain)]
      xTrain_feat <- as.matrix(features[names(xTrain_obs), ])
      if (ncol(xTrain_feat) == 1 ) {
        xTrain_feat <- t(xTrain_feat)
      }
      rownames(xTrain_feat) <- names(xTrain_obs)
      xTrain_feat[is.na(xTrain_feat)] <- 0  ## replace NA by 0
      
      # train a SGL model
      trainOut <- train_SGL(group_index, train_feat, train_obs, xTrain_feat, xTrain_obs, alpha)
      str(trainOut)                       
      
      # extract test set
      test_obs <- res[as.character(crossVal[[i]]$test)]
      test_feat <- as.matrix(features[names(test_obs), ])
      rownames(test_feat) <- names(test_obs)
      test_feat[is.na(test_feat)] <- 0  ## replace NA by 0
      
      # predict with test set
      test_pred <- predictSGL(trainOut$model, test_feat, trainOut$lambda_num ) 
      
      # compute pearson correlation and RMSE
      library(Hmisc) ; library(hydroGOF)
      x <- matrix(c(test_obs, test_pred),nrow=length(test_obs),ncol=2) 
      x <- na.omit(x)
      if( var(x[,1]) ==0) { x[1,1] <- x[1,1] + 1E-100} 
      if( var(x[,2]) ==0) { x[1,2] <- x[1,2] + 1E-100}
      
      # different metrics
      if(length(x[,1]) < 3) {
        pcorr <- 0 # pearson correlation
        pval <- 1 # p value
        result <- matrix(c(pcorr, pval, length(test_obs)),nrow=1,ncol=3) 
        colnames(result) <- c("pcorr", "pval", "Nobs")
        pearson_table <- rbind (pearson_table, result)
      } else {
        
        temp <- cor.test(x[,1], x[,2])
        pcorr <- temp$estimate # pearson correlation
        pval <- temp$p.value # p value
        result <- matrix(c(pcorr, pval, length(test_obs)),nrow=1,ncol=3) 
        colnames(result) <- c("pcorr", "pval", "Nobs")
        pearson_table <- rbind (pearson_table, result)
        
        RMSE_table <- c(RMSE_table , rmse(x[,1], x[,2], na.rm=TRUE) )
     }
    } 
    
  } ## end of iteration  
  result_table <- list(pearson_table, RMSE_table) ; names(result_table) <- c("pcorr","rmse")
  return(result_table)
}

############################################# DEFINE PATHWAY GROUP #############################################

define_group <- function(mat , features) {
 
  mat <- model[rownames(model) %in% colnames(features) , ]
  # mat <- mat[-which(rowSums(mat) == 0) , ] 
  for(i in 1:length(mat[1,])) {  mat[which(abs(mat[,i]) > 0), i] <- i  }
  feature_common <- c()
  for(i in 1:length(mat[,1])) {   
    if( length(which(abs(mat[i, ]) > 0))  > 1 ) {
      feature_common <- c(feature_common , i)
    }
  }
  if(length(feature_common) >= 1) {
  feature_common_names <- rownames(mat)[feature_common]
  mat <- mat[-feature_common, ]
  common_index <- rep(length(colnames(mat))+1, length(feature_common_names))
  common_features <- features[ ,feature_common_names]
  group_ID <- rowSums(mat)
  # OTHER_ID <- colnames(features)[- which(colnames(features) %in% names(group_ID) ) ]
  group_index <- as.numeric(group_ID)
  group_features <- features[ ,names(group_ID)]
  
  index <- c(group_index, common_index)
  features <- cbind(group_features, common_features)
  } else {
  group_ID <- rowSums(mat)
  # OTHER_ID <- colnames(features)[- which(colnames(features) %in% names(group_ID) ) ]
  group_index <- as.numeric(group_ID)
  group_features <- features[ ,names(group_ID)]
  index <- group_index
  features <- group_features
  }
  
  L <- list(index, features)
  return(L)
}
