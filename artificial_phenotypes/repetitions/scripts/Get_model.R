### Gamma tuning
gammaTuning <- function(trainIx, data, target, gamma = 1, ntree = 500){
    res <- lapply(trainIx, wf, data=data, target=target, gamma=gamma, ntree=ntree)
    return(res)
}


### Get the model
wf <- function(ix, data, target, ntree = 500, gamma = 1){
    set.seed(ix[1])
    res <- list()
    # feature selection
    message('Feature selection')
    RF <- RRF(data[ix,], flagReg=0, as.factor(target[ix]))
    imp <-RF$importance[,"MeanDecreaseGini"]
    impRF <- (imp - min(imp))/(max(imp) - min(imp))# normalization
    coefReg <- (1-gamma) + gamma*impRF
    GRRF <- RRF(data[ix,], as.factor(target[ix]), flagReg=1, coefReg=coefReg)
    
    # select data
    message('Subset data')
    to_keep <- unique(c(colnames(data)[GRRF$feaSet], str_subset(colnames(data), pattern='group')))
    X_fs <- select(data, all_of(to_keep))
    
    # RF
    message('RF')
    rf_fs <- randomForest(x = X_fs[ix,], y = target[ix], ntree = ntree)
    pred <- predict(object = rf_fs, newdata = X_fs[-ix, ])
    tmp <- confusionMatrix(data = pred, reference = target[-ix])
    
    return(tmp$overall)
}