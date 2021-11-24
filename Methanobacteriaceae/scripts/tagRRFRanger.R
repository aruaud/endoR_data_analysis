
# feature selection with taxa-aware gRRF and model with ranger
tagRRFRanger <- function(ix, data, meta, target, families, genera, species, gamma = 1, k = 0.5, num.trees = 500){
    
    set.seed(ix[1])
    res <- list()
    data <- mutate_if(data, is.character, as.factor)
    # feature selection
    message('Feature selection')
    RF <- RRF(data[ix,], flagReg=0, as.factor(target[ix]))
    regterm <- data.frame(Feature = names(RF$importance[,"MeanDecreaseGini"])
                     , imp = RF$importance[,"MeanDecreaseGini"])
    # normalization across all features
    regterm$imp_norm <- (regterm$imp - min(regterm$imp))/(max(regterm$imp) - min(regterm$imp))
    # normalization per branch
    regterm$imp_tax <- NA
    regterm$mb <- NA
    for (x in regterm$Feature){
        tmp <- unique(c(families[[x]], genera[[x]], species[[x]]))
        ib <- which(regterm$Feature %in% tmp)
        if (length(ib) == 0) next
        mb <- max(regterm$imp[ib])

        ib <- which(regterm$Feature == x)
        regterm$mb[ib] <- mb
        regterm$imp_tax[ib] <- regterm$imp[ib] / mb
    }
    regterm$imp_tax <- ifelse(is.na(regterm$imp_tax), regterm$imp_norm, regterm$imp_tax)
    # regularization term
    regterm$lambda <- (1-gamma) + gamma*regterm$imp_norm^(1-k)*regterm$imp_tax^k
    GRRF <- RRF(data[ix,], as.factor(target[ix]), flagReg=1, coefReg=regterm$lambda)
    
    # select data
    message('Subset data')
    to_keep <- c(meta, colnames(data)[GRRF$feaSet])
    res$confirmed <- to_keep
    X_fs <- select(data, all_of(unique(to_keep)))
    
    # RF
    message('RF')
    # get the class weights for the weighted forest
    class_weights <- round(sum(target == levels(target)[1])/length(target), digits = 2)
    class_weights <- c(1-class_weights, class_weights)
    
    set.seed(ix[1])
    rf_fs <- ranger(x = X_fs[ix,], y = target[ix], class.weights = class_weights
            , num.trees=num.trees, importance = 'impurity')
    res$rf_model <- rf_fs

    if (length(ix) != nrow(data)){
        pred <- predict(rf_fs, data = X_fs[-ix, ])
        tmp <- confusionMatrix(data = pred$predictions, reference = target[-ix])
        res$rf_performance <- tmp$overall
    }

    return(res)
}


# Ranger results 
getRangerCV <- function(res, ks=seq(0,1, by = 0.1), gammas=seq(0,1, by = 0.1)){
    all <- list()
    for (i in 1:length(res)){
        all[[i]] <- list()
        for (j in 1:length(res[[i]])){
            tmp <- t(sapply(res[[i]][[j]], function(x){return(x$rf_performance)}))
            tmp <- data.frame(meanAcc = mean(tmp[,1]), sdAcc = sd(tmp[,1])
                              , meanK = mean(tmp[,2]), sdK = sd(tmp[,2]), ks = ks[j])
            all[[i]][[j]] <- tmp
        }
        all[[i]] <- do.call(rbind, all[[i]])
        all[[i]]$gamma <- gammas[i]
    }

    all <- do.call(rbind, all)

    return(all %>% arrange(desc(meanK), desc(meanAcc)) )

}
