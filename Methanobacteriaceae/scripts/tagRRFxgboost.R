tagRRFxgboost <- function(ix, data, meta, target, case_weights
    , families, genera, species, gamma = 1, k = 0.5
    , xgboost_param = list(c('nrounds' = 100, 'max_depth' = 6))){
     
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
    
    # XGBoost
    message('XGBoost')
    tuned_xgb <- list()
    for (i in 1:length(xgboost_param)){
        set.seed(ix[1])
        xgb_fs <- xgboost(data = as.matrix(X_fs[ix,]), label = as.character(target[ix])
                        , weight = case_weights[ix]
                        , nrounds = xgboost_param[[i]]['nrounds'], max_depth = xgboost_param[[i]]['max_depth']
                        , objective = 'binary:hinge')
        pred <- predict(object = xgb_fs, newdata = as.matrix(X_fs[-ix,]))
        tmp <- confusionMatrix(data = as.factor(pred), reference = as.factor(target[-ix]))
        tuned_xgb[[i]] <- tmp$overall
    }
    names(tuned_xgb) <- sapply(xgboost_param, paste, collapse = 'nrounds_maxdepth')
    res$tuned_xgb <- tuned_xgb
    
    return(res)
}


# gamma tuning in xgboost
gammaTuning <- function(trainIx, data, meta, target, families, genera, species, case_weights
                        , fs_param = c('gamma' = 1, 'k' = 1)
                        , xgboost_param = list(c('nrounds' = 100, 'max_depth' = 6)) ){
    
    res <- lapply(trainIx, tagRRFxgboost
                  , data=data, meta=meta, target=target, case_weights=case_weights
                  , families = families, genera = genera, species = species
                  , gamma=fs_param['gamma'], k = fs_param['k'], xgboost_param=xgboost_param)
    
    return(res)
}




# xgboost results 
getXGBoostCV <- function(res, ks=seq(0,1, by = 0.1), gammas=seq(0,1, by = 0.1)){
    all <- list()
    for (i in 1:length(res)){
        all[[i]] <- list()
        for (j in 1:length(res[[i]])){
            tmp <- t(sapply(res[[i]][[j]], function(x){return(x$xgb_performance)}))
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
