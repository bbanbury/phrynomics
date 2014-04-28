CalculateMissingData <- function(data){
  splitdata <- SplitSNP(data)
  MissingSites <- apply(splitdata, 2, function(x) length(which(x == "N")))
  PercentMissing <- MissingSites/dim(splitdata)[1]
  return(PercentMissing)
}