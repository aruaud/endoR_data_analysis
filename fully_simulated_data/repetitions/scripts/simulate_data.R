getSimulations <- function(nr = 5000, v_tosamp, prandom = 0.2, seedOri=1){
    simu <- data.table(tmp = 1:nr)
    
    # get the vectors
    for (j in 1:12){
        set.seed(seedOri*j + seedOri)
        set(simu, NULL, paste0('V',j), value = sample(v_tosamp, size = nr, replace = TRUE))
    }
    simu <- simu[,tmp:=NULL]
    
    # grpoup row numbers
    i_a <- 1:(0.25*nr)
    i_b <- (0.25*nr+1):(0.5*nr)
    i_c <- (0.5*nr+1):(0.75*nr)
    i_d <- (0.75*nr+1):nr
    
    # get the target
    set(simu, i_a, 'target', sign(simu[[1]][ i_a ]*simu[[2]][ i_a ]) )
    set(simu, i_a, 'group', 'a')

    set(simu, i_b, 'target', sign(simu[[3]][ i_b ]) )
    set(simu, i_b, 'group', 'b')

    set(simu, i_c, 'target', ifelse(simu[[4]][i_c]>0&simu[[5]][i_c]>0, 1, -1))
    set(simu, i_c, 'group', 'c')

    set(simu, i_d, 'target', ifelse(simu[[6]][i_d]>0&simu[[7]][i_d]<0, 1, -1))
    set(simu, i_d, 'group', 'd')
    
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