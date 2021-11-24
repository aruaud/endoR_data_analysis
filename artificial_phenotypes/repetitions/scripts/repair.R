repairEdges <- function(f, X, rules,related_taxa, K){
	ori <- qread(f)
	seedOri <- as.numeric(str_extract(f, pattern = '(?<=simu)[:digit:]+(?=\\_)'))
	target <- ori$target

	## Randomize samples
	set.seed(seedOri)
	X <- X[sample(1:nrow(X)),]

	#### Get the exec
	set.seed(seedOri)
	# first sampling
	rownames(rules) <- 1:nrow(rules)
	i_n <- sample(1:nrow(rules), 1)
	# second sampling: remove the rules with var already used
	tmp <- rownames(subset(rules, !(v1 %in% c(rules$v1[i_n], rules$v2[i_n])) & !(v2 %in% c(rules$v1[i_n], rules$v2[i_n]))  ))
	i_n <- c(i_n, sample(tmp, 1))
	# last sampling: remove the rules with var already used
	tmp <- rownames(subset(rules, !(v1 %in% c(rules$v1[i_n], rules$v2[i_n])) & !(v1 %in% c(rules$v1[i_n], rules$v2[i_n]))  ))
	i_n <- c(i_n, sample(tmp, 1))	   
	exec <- rules[as.numeric(i_n), 'condition']

	# transform to dummy
	X$group <- as.factor(ori$group)
	dummies <- dummyVars(~ ., data = X )
	dummies <- as.data.table(predict(dummies, newdata = X ))
	colnames(dummies) <- colnames(dummies) %>% str_replace_all(pattern = '\\.', replacement = '')


	groups <- c('a', 'b', 'c', 'd')


	### Get the full rules and all
	exec <- unlist(lapply(1:3,getAllRules, exec = exec, groups= groups, colN=colnames(dummies)))
	exec <- data.table(condition = exec, n = 1) 
	exec <- changeRulesDummies(rules = exec, dummy_var = 'group', data = dummies
	                           , target =target, classPos = '1', in_parallel = FALSE)
	exec <- discretizeRules(rules = exec, data = dummies, target=target, K = K, classPos='1'
	                       , in_parallel = TRUE, n_cores = 2)
	data_ctg <- exec$data_ctg
	exec <- exec$rules_ctg

	mod <- model2DE(data = data_ctg, target = target, classPos = '1'
	                    , exec = exec, filter = FALSE, prune = TRUE, maxDecay = 0.0001, in_parallel = FALSE)
	true_edges <- select(mod$edges, c('x', 'y'))
	true_edges <- asplit(as.matrix(true_edges), MARGIN = 1)        
	    
	#### Get the edges lists
	# this looks ugly but not everything can look nice, it is functional at least
	# the expanded one
	related_taxa$groupa <- 'groupa'
	related_taxa$groupb <- 'groupb'
	related_taxa$groupc <- 'groupc'
	    
	expanded_edges <- list()
	for (i in 1:length(true_edges)){
	    suf <- true_edges[[i]] %>% str_extract_all(pattern = '\\_{2}.*')
	    tmp <- true_edges[[i]] %>% str_replace(pattern = '\\_{2}.*', replacement = '')

	    tmp <- true_edges[[i]] %>% str_replace(pattern = '\\_{2}.*', replacement = '')
	    tmp <- related_taxa[tmp]
	    tmp <- expand.grid(paste0(tmp[[1]], suf[[1]]), paste0(tmp[[2]], suf[[2]]))  
	    tmp <- asplit(tmp, MARGIN=1)
	    
	    expanded_edges <- append(expanded_edges, tmp)
	}

	ori$true_edgesK2 <- true_edges
	ori$expanded_edgesK2 <- expanded_edges
	ori$exec <- exec
	ori$rules <- rules[as.numeric(i_n), ]

	return(ori)

	#message('Save the object')
    #qsave(x = ori, file = f)

    #return('success!')
}