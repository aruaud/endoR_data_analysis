# From Hadley Whickham, https://r.789695.n4.nabble.com/Suppressing-output-e-g-from-cat-td859876.html
quiet <- function(x) {
  sink(tempfile())
  on.exit(sink())
  invisible(force(x))
} 

# Calculate kappa from a confusion matrix
getKappa <- function(x){
    p0 <- (x[1,1] + x[2,2]) / sum(x)
    py <- (x[1,1] + x[1,2])*(x[1,1] + x[2,1]) / sum(x)^2
    pn <- (x[2,1] + x[2,2])*(x[1,2] + x[2,2]) / sum(x)^2
    pe <- py+pn
    
    k <- 1 - (1-p0)/(1-pe)
    return(k)
}