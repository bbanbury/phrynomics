RemoveNonBinary <- function(data){
  splitdata <- SplitSNP(data)
  binaryVector <- apply(splitdata, 2, IsBinary)
  data2 <- cSNP(splitdata, KeepVector=binaryVector)
  return(data2)
}
