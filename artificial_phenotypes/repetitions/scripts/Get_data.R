### Get data
getResampleData <- function(X,rules,related_taxa,K, prandom = 0.2, seedOri=1){
# rules: the corresponding rules
# seedOri: the seed
# related_taxa: list of related taxa
#for a family: all genera and species, 
#for a genus: the family and all species,
#for a species: the family and genus, and all species in the same genus
    
## Randomize samples
set.seed(seedOri)
X <- X[sample(1:nrow(X)),]

#### Get the target
set.seed(seedOri)
# first sampling
rownames(rules) <- 1:nrow(rules)
i_n <- sample(1:nrow(rules), 1)
# second sampling: remove the rules with var already used
tmp <- rownames(subset(rules, !(v1 %in% c(rules$v1[i_n], rules$v2[i_n])) & !(v2 %in% c(rules$v1[i_n], rules$v2[i_n]))  ))
i_n <- c(i_n, sample(tmp, 1))
# last sampling: remove the rules with var already used
tmp <- rownames(subset(rules, !(v1 %in% c(rules$v1[i_n], rules$v2[i_n])) & !(v2 %in% c(rules$v1[i_n], rules$v2[i_n]))  ))
i_n <- c(i_n, sample(tmp, 1))
    
exec <- rules[as.numeric(i_n), 'condition']

target <- data.frame(Sample = 1:nrow(X))
target$groupa <- eval(parse(text = exec[1] ))
target$groupb <- eval(parse(text =  exec[2] ))
target$groupc <- eval(parse(text = exec[3] ))

# make the groups
set.seed(seedOri)
nr <- nrow(X)
ng <- floor(nr/3)
X$group <- c(rep('a',ng), rep('b', ng), rep('c', nr-2*ng))
target_c <- c(target$groupa[X$group == 'a'], target$groupb[X$group == 'b'], target$groupc[X$group == 'c'])

# add noise
groups <- c('a', 'b', 'c', 'd')
set.seed(seedOri)
brnounou <- rbinom(n = length(X$group), size = 1,prob = prandom)
for (i in 1:length(brnounou)){
    if (brnounou[i] == 1){
        set.seed(i)
        X$group[i] <- sample(groups[groups != X$group[i]], 1)
    }
}

# format
target_c <- as.factor(ifelse(target_c == TRUE, '1', '-1'))
X$group <- as.factor(X$group)
# transform to dummy
dummies <- dummyVars(~ ., data = X )
dummies <- as.data.table(predict(dummies, newdata = X ))
colnames(dummies) <- colnames(dummies) %>% str_replace_all(pattern = '\\.', replacement = '')


### Get the full rules and all
exec <- unlist(lapply(1:3,getAllRules, exec = exec, groups= groups, colN=colnames(dummies)))
exec <- data.table(condition = exec, n = 1) 
exec <- changeDecisionsDummies(rules = exec, dummy_var = 'group', data = dummies
                           , target =target_c, classPos = '1', in_parallel = FALSE)
exec <- discretizeDecisions(rules = exec, data = dummies, target=target_c
                  , K = K, classPos='1'
                  , in_parallel = TRUE, n_cores = 2)
data_ctg <- exec$data_ctg
exec <- exec$rules_ctg
mod <- model2DE(data = data_ctg, target = target_c, classPos = '1'
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

return(list('x' = dummies, 'y' = target_c
     , 'group' = X$group, 'true_edges' = true_edges, 'expanded_edges' = expanded_edges
     #, 'ground_truth' = mod
     #, 'exec' = exec
           ))

}



getAllRules <- function(i, exec, groups, colN){
    gp <- paste0('group', groups[i])
    gpR <- paste0('X[,',which(colN == gp), ']>0.5')
    subR <- unlist(strsplit(str_replace_all(exec[i], pattern = '>', replacement = '<='), split = ' & '))
    
    res <- c(paste(gpR, exec[i], sep = ' & '), paste(gpR, subR, sep = ' & '))
    return(res)
}