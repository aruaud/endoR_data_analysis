### Functions for the run of endoR on different bootstraps for a same replicate.

## Wrapper function for the resampling: get the sample weight for bootstrapping and save the object at the end
wrappy <- function(seedN = 0, rf_fs, target, data
                   , path = '/ebio/abt3_projects/temp_data/aruaud/MtgSimu50/p005_B100-10same_2/'
                   , times = 100, p=0.5
                   , discretize = TRUE, K = 2, prune = TRUE
                   ,in_parallel = TRUE, n_cores = 10){
    
    pPos <- round(sum(target == '1')/length(target), digits = 1)
    sample_weight <- c(pPos, 1-pPos)
    
    # network
    message('Run endoR')
    resampled <- 
        model2DE_resampling(model = rf_fs, model_type = 'rf', dummy_var = 'group'
                     ,times = times, p = p, sample_weight = sample_weight
                     ,ntree = 'all' 
                     ,data = data, target = target, classPos = '1'
                     ,discretize = discretize, K = K
                     ,prune = prune, seed = seedN
                     ,in_parallel = in_parallel, n_cores = n_cores
                    ) 
    
    message('Save the object')
    fname <- paste0(path, 'simu',seedN, '.qs')
    qsave(x = resampled$resamp, file = fname)
    return('success!')
}


### Get one precision/recall curve 
getSinglePR_singleB <- function(pc, expanded_edges, all){
    # prepare edges from endoR
    edges <- pc$edges_agg %>% select(x, y, condition, importance, influence, d.x, d.y) %>% 
        subset(condition %in% pc$rules_summary$condition) %>%
        group_by(x, y, condition) %>% 
        summarise(importance= mean(importance), influence = mean(influence*(d.x+d.y)/2)) %>%
        ungroup %>% 
        left_join(select(pc$rules_summary,condition, inN, imp, n), by = 'condition') %>% 
        select(-condition) 
    
    
    # get the PR
    edi <- lapply(unique(edges$inN)
            , function(i, edges, true_edges,expanded_edges, all){
                suppressMessages(edges <- edges %>% subset(inN >= i) %>% group_by(x, y) %>% 
                                 summarise( inN = max(inN), importance = sum(importance*imp*n)
                                           , influence = sum(influence*imp*n)/sum(imp*n) ) %>% ungroup)
                metricsNet(edges = edges,expanded_edges=expanded_edges, all=all)
              }, edges=edges,expanded_edges=expanded_edges, all=all)
    edi <- as.data.frame(do.call(rbind, edi))
    edi$inN <- unique(edges$inN)
    edi <- add_row(edi, n_edges=0, tp=0,fp=0,tn=0,fn=0,inN=0)
    
    return(edi %>% arrange(tp))
}


### Get all precision/recall curves from the same replicate on different bootstraps
getPR_singleB <- function(res, alpha, pi_thr, related_taxa, all){
    # get the true edges
    te <- lapply(res$true_edges, str_replace_all, pattern = '\\_{2}.*', replacement = '')
    te <- unique(lapply(te, sort))
    expanded_edges <- list()
    for (i in 1:length(te)){
        tmp <- te[[i]] %>% str_replace(pattern = '\\_{2}.*', replacement = '')
        tmp <- related_taxa[tmp]
        tmp <- expand.grid(tmp[[1]], tmp[[2]])  
        tmp <- asplit(tmp, MARGIN=1)

        expanded_edges[[i]] <- sapply(lapply(tmp, sort), paste, collapse = ' - ')
        names(expanded_edges)[i] <- paste(sort(te[[i]]), collapse = ' - ')
    }
    
    # other parameters
    pc <- stabilitySelection(res$resampled, alpha_error = alpha, pi_thr = pi_thr)
    edi <- getSinglePR_singleB(pc, expanded_edges, all = all)
    edi$B <- length(res$resamp)
    edi$K <- length(expanded_edges)
    
    return(edi)
}


### Interpolate curves
getInterpo <- function(i, raw){
    raw <- subset(raw, rep == i)
    res <- list()
    
    N <- unique(raw$Ncol)
    N <- (N^2-N)/2
    K <- unique(raw$K)
    
    j <- 1
    for (a in unique(raw$B)){
        tmp <- approx(raw$fp[raw$B == a]
                       , raw$tp[raw$B == a]
                       , xout = seq(1, 144, length.out = 1000)
                       , ties = max)
        tmp <- as.data.frame(do.call(cbind, tmp))
        colnames(tmp) <- c('fp', 'tp')
        
        # add 0 and max values to the extrapolations
        #mifp <- min(tmp$fp[!is.na(tmp$tp)])
        mafp <- max(tmp$fp[!is.na(tmp$tp)])
        #tmp$tp[tmp$fp<mifp] <- min(tmp$tp)
        tmp$tp[tmp$fp>mafp] <- max(tmp$tp, na.rm = TRUE) #K
        
        tmp <- tmp[complete.cases(tmp),]
        #tmp <- tmp %>% add_row(tp = 0, fp = 0)
        tmp$B <- a

        res[[i]] <- tmp
        j <- j+1
    }
    res <- do.call(rbind, res)
    res$rep <- i
    
    return(res)
}