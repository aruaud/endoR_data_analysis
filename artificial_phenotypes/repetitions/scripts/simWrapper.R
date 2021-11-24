# Wrapper
simWrapper <- function(
    data,rules, related_taxa
    , seedOri, prandom = 0.2
    , gammas = seq(0,1, by = 0.1), ntree=500
    , times = 5, p = 0.75
    , discretize = TRUE, K = 3,prune = TRUE
    , in_parallel = TRUE, n_cores = 5
    , path = '/ebio/abt3_projects/temp_data/aruaud/'
){
    res <- list()
    #on.exit(return(res))
    
    #simulations
    simu <- getResampleData(X = data, rules=rules, related_taxa=related_taxa, K = K
                            , prandom = prandom, seedOri=seedOri)
    X <- simu$x
    target <- simu$y
    res$target <- target
    res$expanded_edges <- simu$expanded_edges
    res$true_edges <- simu$true_edges
    res$group <- simu$group


    message('Tuning')
    # gamma tuning
    set.seed(seedOri)
    trainIx <- createDataPartition(y = target, times = 10, p = .7, list = TRUE)
    gam_tun <- lapply(gammas, gammaTuning, data = X, target = target, trainIx = trainIx, ntree = ntree)
    gam_tun <- gam_tun %>% lapply(do.call, what = rbind) %>% 
            sapply(function(x){return(c('meanAcc' = mean(x[,1]*100), 'sdAcc' = sd(x[,1]*100)
                                        , 'meanK' = mean(x[,2]), 'sdK' = sd(x[,2])))})
    gam_tun <- t(gam_tun)
    rownames(gam_tun) <- gammas
    gam_tun <- as.data.frame(gam_tun) %>% arrange(desc(meanK), desc(meanAcc))
    gamma <- as.numeric(rownames(gam_tun)[1])
    res$gamma_all <- gam_tun
    res$gamma <- gamma
    
    message('Fit final model')
    #model
    set.seed(seedOri)
    RF <- RRF(X, flagReg=0, target)
    imp <- RF$importance[,"MeanDecreaseGini"]
    impRF <- (imp - min(imp))/(max(imp) - min(imp))# normalization
    coefReg <- (1-gamma) + gamma*impRF
    GRRF <- RRF(X, target, flagReg=1, coefReg=coefReg)
    # select data
    to_keep <- unique(c(colnames(X)[GRRF$feaSet], str_subset(colnames(X), pattern='group')))
    X_fs <- select(X, all_of(to_keep))
    res$data <- X_fs
    # RF
    set.seed(seedOri)
    rf_fs <- randomForest(x = X_fs, y = target, ntree = ntree)
    
    # model metrics
    res$rf_accuracy <- (rf_fs$confusion[1,1] + rf_fs$confusion[2,2])/nrow(X_fs)
    res$rf <- rf_fs

    pPos <- round(sum(simu$y == '1')/length(simu$y), digits = 1)
    sample_weight <- c(pPos, 1-pPos) 


    # network
    message('Run endoR')
    resampled <- 
        model2DE_resampling(model = rf_fs, model_type = 'rf', dummy_var = 'group'
                     ,times = times, p = p, sample_weight = sample_weight
                     ,ntree = 'all' 
                     ,data = X_fs, target = target, classPos = '1'
                     ,discretize = discretize, K = K
                     ,prune = prune
                     ,in_parallel = in_parallel, n_cores = n_cores
                    ) 
    res$resampled <- resampled$resamp
    
    # get the column names ... don't ask why it is so complicated..
    message('Get the column names')
    tmp <- discretizeData(X_fs, K = K, return_split = FALSE) 
    notNum <- which(!sapply(tmp, is.numeric))
    for (j in notNum) set(tmp,i = NULL ,j, paste0('__', tmp[[j]]))
    data_ctg <- dummyVars( ~ ., data = tmp)
    data_ctg <- as.data.table(predict(data_ctg, newdata = tmp))
    for (j in notNum) set(data_ctg,i = NULL ,j, str_replace(data_ctg[[j]], pattern = '^__', replacement = ''))
    colnames(data_ctg) <- gsub(x = colnames(data_ctg), pattern = "`", replacement="")
    colN <- colnames(data_ctg)
    res$colN <- colN
    
    #network accuracy 
    #message('Get network accuracy and so on')
    #nacc <- list()
    #i <- 1
    #for (m in minN){
    #    tmp <- sapply(alphas, getNetMeasures, rules = resampled$resamp
    #               , minN = m, true_edges=simu$true_edges, expanded_edges=simu$expanded_edges, colN = colN)
    #    tmp <- as.data.frame(t(tmp))
    #    tmp$alpha <- alphas
    #    tmp$p_thr <- m
    #    nacc[[i]] <- tmp
    #    i <- i+1
    #}
    #nacc <- do.call(rbind, nacc)
    #res$network_acc <- nacc
    
    message('Save the object')
    prandom <- str_replace(as.character(prandom), pattern = '\\.', replacement = '_')
    fname <- paste0(path, 'simu',seedOri,'_p', prandom, '.qs')
    qsave(x = res, file = fname)
    
    
    return('success!')
}
