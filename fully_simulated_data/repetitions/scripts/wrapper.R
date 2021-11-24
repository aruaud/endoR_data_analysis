simWrapper <- function(
    seedOri, v_tosamp, nr=5000, prandom = 0.2
    , times = 5, p = 0.75, ntree = 500
    , discretize = TRUE, K = 3,prune = TRUE
    , in_parallel = TRUE, n_cores = 5
    , path = '/ebio/abt3_projects/temp_data/aruaud/'
){
    res <- list()
    
    # simulations
    simu <- getSimulations(nr = nr, v_tosamp=v_tosamp, prandom=prandom, seedOri=seedOri)
    res$data <- simu

    # model
    set.seed(0)
    rf <- randomForest(y = simu$y, x = simu$x, ntree = ntree)
    
    # model metrics
    res$rf_accuracy <- (rf$confusion[1,1] + rf$confusion[2,2])/nr
    tmp <- as.data.frame(rf$importance)
    tmp$Feature <- factor(rownames(tmp), levels = rownames(tmp)[order(tmp$MeanDecreaseGini)])
    res$FeatureGini <- arrange(tmp, desc(MeanDecreaseGini))
    
    pPos <- round(sum(simu$y == '1')/length(simu$y), digits = 1)
    sample_weight <- c(pPos, 1-pPos) 
    
    # network
    resampled <- 
        model2DE_resampling(model = rf, model_type = 'rf', dummy_var = 'group'
                     ,times = times, p = p, sample_weight = sample_weight
                     ,ntree = 'all' 
                     ,data = simu$x, target = simu$y, classPos = '1'
                     ,discretize = discretize, K = K
                     ,prune = prune, filter = FALSE
                     ,in_parallel = in_parallel, n_cores = n_cores
                    ) 
    res$resamp <- resampled$resamp

    fname <- paste0(path, 'simu',seedOri,'_p', floor(prandom*10),'_n', nr, '.qs')
    qsave(x = res, file = fname)
    
    return('rf_accuracy'= res$rf_accuracy)
}