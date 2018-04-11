####################################################################### randomForest ##############################################################################
###################################################################################################################################################################
library(randomForest)
library(hydroGOF)
 
train_RF <- function(res, train_percent=train_percent) {
  RES_features <- cbind(res, features) ; colnames(RES_features)[1] <- "Res"
  RES_features <- RES_features[complete.cases(RES_features), ]
  train <- sample(1:nrow(RES_features), round(nrow(RES_features) * train_percent))
  RES_features.obj <- randomForest(Res ~ ., data = RES_features[train, ],ntree=100)
  RES_features.pred <- predict(RES_features.obj, RES_features[-train , ])
  
  corr <- cor.test(RES_features[-train , 1],RES_features.pred)
  pcorr <- corr$estimate
  rmse <- rmse(RES_features[-train , 1],RES_features.pred )
  result <- c(pcorr,rmse) ; names(result) <- c("pcorr","rmse")
  return(result)
}

run_RF <- function(index, RES,features, train_percent=train_percent, iteration) {
  RES <- as.matrix(RES)
  features <- as.matrix(features)
  
  res <- RES[index, ]
  result_table <- c()
  for(i in 1:iteration) { result_table <- rbind(result_table, train_RF(res,train_percent=train_percent) ) }
  
  result <- colMeans(result_table)
  result_sd <- apply(result_table,2,sd)
  
  result_table <- c( result[1],result_sd[1] , result[2],result_sd[2])
  names(result_table) <- c("pcorr", "pcorr_sd", "rmse", "rmse_sd")
  return(result_table)
}

SAVE_MODEL_RF <- function(res, train_percent=train_percent ) {
  
  RES_features <- cbind(res, features) ; colnames(RES_features)[1] <- "Res"
  RES_features <- RES_features[complete.cases(RES_features), ]
  train <- sample(1:nrow(RES_features), round(nrow(RES_features) * train_percent))
  RES_features.obj <- randomForest(Res ~ ., data = RES_features[train, ], ntree = 100)
  
  return(RES_features.obj)
  
}


# ####################################################################### randomForestSRC ###########################################################################
# ###################################################################################################################################################################
# library(randomForestSRC)
# library(hydroGOF)
# 
# train_RF_SRC <- function(res, train_percent=train_percent) {
#   RES_features <- cbind(res, features) ; colnames(RES_features)[1] <- "Res"
#   RES_features <- RES_features[complete.cases(RES_features), ]
#   train <- sample(1:nrow(RES_features), round(nrow(RES_features) * train_percent))
#   RES_features.obj <- rfsrc(Res ~ ., data = RES_features[train, ])
#   RES_features.pred <- predict(RES_features.obj, RES_features[-train , ])
#   
#   corr <- cor.test(RES_features.pred$yvar,RES_features.pred$predicted)
#   pcorr <- corr$estimate
#   rmse <- rmse(as.numeric(RES_features.pred$yvar),as.numeric(RES_features.pred$predicted) )
#   result <- c(pcorr,rmse) ; names(result) <- c("pcorr","rmse")
#   return(result)
# }
# 
# run_RF_SRC <- function(index, RES,features, train_percent=train_percent, iteration) {
#   res <- RES[i, ]   
#   result_table <- c()
#   for(i in 1:iteration) { result_table <- rbind(result_table, train_RF_SRC(res,train_percent=train_percent) ) }
#   
#   result <- colMeans(result_table)
#   result_sd <- apply(result_table,2,sd)
#   
#   result_table <- c( result[1],result_sd[1] , result[2],result_sd[2])
#   names(result_table) <- c("pcorr", "pcorr_sd", "rmse", "rmse_sd")
#   return(result_table)
# }
# 
# SAVE_MODEL_RF_SRC <- function(res, train_percent=train_percent ) {
#   
#   RES_features <- cbind(res, features) ; colnames(RES_features)[1] <- "Res"
#   RES_features <- RES_features[complete.cases(RES_features), ]
#   train <- sample(1:nrow(RES_features), round(nrow(RES_features) * train_percent))
#   RES_features.obj <- rfsrc(Res ~ ., data = RES_features[train, ])
#   
#   return(RES_features.obj)
# 
# }
