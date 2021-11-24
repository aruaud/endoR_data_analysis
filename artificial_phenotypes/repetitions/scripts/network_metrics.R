### Get the network accuracies
metricsNet <- function(edges, expanded_edges, all){
    if (nrow(edges) == 0){
        return(c('n_edges' = 0, 'tp' = NA, 'fp' = NA, 'tn' = NA, 'fn' = NA))
    } else { n_edges <- nrow(edges) }
    
    pred <- unique(select(edges, c('x', 'y'))) 
    pred <- asplit(as.matrix(pred), MARGIN = 1)
    pred <- lapply(pred, str_replace, pattern = '\\_{2}.*', replacement = '')
    pred <- unique( lapply(pred, sort) )
    pred_edges <- sapply(pred, function(x){paste(x, collapse = ' - ')})    
    
    
    # those that should not be but are = in pred but not truth
    fp <- sum(sapply(pred_edges, function(x){!(x %in% unlist(expanded_edges))}))
    # those that should be and are
    tp <- sum(sapply(expanded_edges, function(x){ifelse(sum(x %in% pred_edges) > 0, 1, 0)}))
    # those that should be but are not = in truth but not in pred
    fn <- sum(sapply(names(expanded_edges), function(x){!(x %in% pred_edges)}))
    #
    tn <- all - ( fp+tp+fn )

    res <- c('n_edges' = n_edges, 'tp' = tp, 'fp' = fp, 'tn' = tn, 'fn' = fn) 
    return(res)

}


getNetMeasures <- function(rules, alpha_error, minN,expanded_edges, all){
    pc <- postCluster(res = rules, alpha_error = alpha_error, minN = minN)
    n_decisions <- nrow(subset(pc$rules_summary, inN >= length(rules)*minN))
    metnet <- round(metricsNet(edges = pc$edges,expanded_edges=expanded_edges, all = all), digits = 4)
    
    return(c('n_decisions' = n_decisions, metnet))
}
