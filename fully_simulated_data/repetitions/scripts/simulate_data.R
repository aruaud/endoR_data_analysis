getSimulations <- function(nr = 5000, v_tosamp, prandom = 0.2, seedOri=1){
    simu <- data.table(tmp = 1:nr)
    
    # get the vectors
    for (j in 1:12){
        set.seed(seedOri*j + seedOri)
        set(simu, NULL, paste0('V',j), value = sample(v_tosamp, size = nr, replace = TRUE))
    }
    simu <- simu[,tmp:=NULL]
    
    # get the target
    set(simu, 1:(0.25*nr), 'target', simu[[1]][ 1:(0.25*nr) ]*simu[[2]][ 1:(0.25*nr) ] )
    set(simu, 1:(0.25*nr), 'group', 'a')

    set(simu, (0.25*nr+1):(0.5*nr), 'target', simu[[3]][ (0.25*nr+1):(0.5*nr) ])
    set(simu, (0.25*nr+1):(0.5*nr), 'group', 'b')

    set(simu, (0.5*nr+1):(0.75*nr), 'target', simu[[4]][(0.5*nr+1):(0.75*nr) ]+simu[[5]][ (0.5*nr+1):(0.75*nr) ])
    set(simu, (0.5*nr+1):(0.75*nr), 'group', 'c')

    set(simu, (0.75*nr+1):nr, 'target', simu[[6]][(0.75*nr+1):nr] - simu[[7]][(0.75*nr+1):nr])
    set(simu, (0.75*nr+1):nr, 'group', 'd')
    
    # randomise
    if (prandom > 0){
        set.seed(seedOri)
        brnounou <- rbinom(n = nr, size = 1,prob = prandom)
        simu$target[brnounou == 1] <- -simu$target[brnounou == 1]
    }
    
    y <- as.factor(sign(simu$target))
    
    # transform to dummy
    dummies <- dummyVars(~ ., data = simu )
    dummies <- as.data.table(predict(dummies, newdata = simu ))
    dummies <- select(dummies, -target)
    
    return( list('x' = dummies, 'y' = y) )

}