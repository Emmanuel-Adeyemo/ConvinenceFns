# For calculating ARVs especially for genesyt exports

arv <- function(dat){
  
  cols <- colnames(dat)
  for(i in 1:ncol(dat)){
    
    x <- cols[i]
    if(is.numeric(dat[,i])){
      dat[sapply(dat, is.numeric)]
      x_name <- paste0(x, "_arv")
      dat[,x_name] <-  100+(sapply(dat[,x], mean) - mean(dat[,x]))/sd(dat[,x])*10
    }
  }
  
 
  return(dat)
}

arv(test)
