metricsNet <- function(edges, true_edges, colN){
    if (nrow(edges) == 0){
        return(c('n_edges' = 0, 'tp' = NA, 'fp' = NA, 'tn' = NA, 'fn' = NA)) #, 'tp_imp' = NA, 'fp_imp' = NA
    } else { n_edges <- nrow(edges) }
    
    res <- numeric()
    # truth: list of vectors of length 2 with the true pairs of associated variables
    # pred: same but as predicted by the method
    all <- asplit(combn(colN, m = 2), MARGIN = 2)
    
    pred <- unique(select(edges, c('x', 'y'))) #, 'importance'
    pred <- asplit(as.matrix(pred), MARGIN = 1)
    pred <- lapply(pred, str_replace, pattern = '\\_{2}.*', replacement = '')
    pred <- unique( lapply(pred, sort) )
    pred_edges <- sapply(pred, function(x){paste(x, collapse = ' - ')})    
    
    true_edges <- lapply(true_edges, sort)
    truth <- sapply(true_edges, function(x){paste(x, collapse = ' - ')})
    
    # those that should not be but are = in pred_edges but not truth
	tp <- length(which(pred_edges %in% truth))
	# those that should be and are
	fp <- length(which(!(pred_edges %in% truth)))
	# those that should be but are not = in truth but not in pred_edges
	fn <- length(which(!(truth %in% pred_edges)))
	#
	tn <- length(all) - ( fp+tp+fn )
    
    # sum of FP importances
    #if(length(fp) > 0) {
    #    fp_imp <- sum(sapply(pred[fp], function(x){as.numeric(x['importance'])}))
    #} else {fp_imp <- 0}
    # sum of TP importances
    #if(length(tp) > 0) {
    #    tp_imp <- sum(sapply(pred[tp], function(x){as.numeric(x['importance'])}))
    #} else {tp_imp <- 0}
    
    res <- c('n_edges' = n_edges, 'tp' = tp, 'fp' = fp, 'tn' = tn, 'fn' = fn) #, 'tp_imp' = tp_imp, 'fp_imp' = fp_imp
    return(res)

    #accuracy:
    acc <- ( length(tp) + tn )/length(all)
    sens <- length(tp)/(length(tp)+fn)
    spec <- tn/(tn+length(fp))
    prec <- length(tp)/(length(tp)+length(fp))
    fdr <- length(fp)/(length(fp)+length(tp))
    # with weighted by the importance (possible only for drawn = positive edges)
    fdr_imp <- fp_imp/(fp_imp+tp_imp)
    prec_imp <- tp_imp/(fp_imp+tp_imp)
    
    
}

getNetMeasures <- function(rules, alpha_error, minN, true_edges, colN){
    pc <- postCluster(res = rules, alpha_error = alpha_error, minN = minN)
    n_decisions <- nrow(subset(pc$rules_summary, inN >= length(rules)*minN))
    metnet <- round(metricsNet(edges = pc$edges, true_edges = true_edges, colN = colN), digits = 4)

    return(c('n_decisions' = n_decisions, metnet))
}
