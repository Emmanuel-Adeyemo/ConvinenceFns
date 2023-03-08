# For calculating ARVs especially for genesyt exports

library(crayon)
library(ggplot2)

#'@param dat - dataset 
#'@param colums_for_arv - column names of interest for calc ARV

do.arv <- function(x, dat){
  
  if(!is.numeric(dat[,x])){
    stop(cat(blue("The specified column name -", x, "- is not numeric.\n")))
  }
  
  out_col_name <- paste0(x, "_arv")
  solve_col <-  round((100+(sapply(dat[,x], mean) - mean(dat[,x]))/sd(dat[,x])*10),0)
  solve_col <- as.data.frame(solve_col)
  colnames(solve_col) <- out_col_name
  
  return(solve_col)
}


out.arv <- function(dat, columns_for_arv){
  
  arv_res <- lapply(columns_for_arv, do.arv, dat)
  
  arv_bind <- bind_cols(arv_res)
  
  dt_arv <- cbind(dat, arv_bind)
}


### example with mpg dataset from the ggplot2 library
d <- as.data.frame(mpg)

sapply(mpg, class)

col_of_interest <- c("displ", "cty", "hwy")

ex <- out.arv(d, col_of_interest)

