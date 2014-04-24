MakePresentAbsent <- function(dataset){
#function removes loci breaks, and changes all data to 1, all missing data to 0
  split <- SplitSNP(dataset)
  split <- split[,-which(split[1,] == " ")]
  split[c(which(split=="-"), which(split=="N"))] <- 0
  split[-c(which(split=="-"), which(split=="N"), which(split=="0"))] <- 1
  return(matrix(as.numeric(unlist(split)),nrow=nrow(split)))
}
