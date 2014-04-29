CalculateMissingData <- function(data, missingChar="N"){
  splitdata <- SplitSNP(data)
  MissingSites <- apply(splitdata, 2, function(x) length(which(x == missingChar)))
  PercentMissing <- MissingSites/dim(splitdata)[1]
  return(PercentMissing)
}
