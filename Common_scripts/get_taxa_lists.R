# Functions to get lists of finer taxa for each family and genus

getFamilies <- function(x, tax_names){
    tmp <- unique(unlist(subset(tax_names, f == x, select = c(f,g))))
    return(tmp[!is.na(tmp)])
}

getGenera <- function(x, tax_names){
    tmp <- unique(unlist(subset(tax_names, g == x, select = c(f,g,s))))
    return(tmp[!is.na(tmp)])
}


# Function to get all other species from the same genus 
getSpecies <- function(x, tax_names){
    tmp <- unique(tax_names$g[tax_names$s == x])
    tmp <- tmp[!is.na(tmp)]
    tmp <- unique(unlist(subset(tax_names, g == tmp, select = c(g,s))))
    return(tmp[!is.na(tmp)])
}