##' Calculate Amount of Missing Data
#' 
#' This function will calculate what percentage of the dataset is missing.  You must give it the missing character value (ex: "N", "-"), and you can not have more than one value in the dataset. 
#' @param Data SNP data in the class "matrix", "data.frame", or "snp"
#' @param missingChar The character that represents missing data
#' @export
#' @return Returns a vector of missing data for each SNP
#' @seealso \link{ReadSNP} \link{ConvertMissingData}
#' @examples
#' data(fakeData)
#' CalculateMissingData(fakeData)

CalculateMissingData <- function(data, missingChar="N"){
  if(class(data) == "snp")
    data <- data$data 
  splitdata <- SplitSNP(data)
  MissingSites <- apply(splitdata, 2, function(x) length(which(x == missingChar)))
  PercentMissing <- MissingSites/dim(splitdata)[1]
  return(PercentMissing)
}
