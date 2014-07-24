##' Get Base Frequencies
#' 
#' This function will calculate base frequencies for any single site. It will not return missing or ambiguous data, so values may not add up to one.  
#' @param site Single SNP site as a vector
#' @export
#' @return Returns a named vector of frequencies.
#' @seealso \link{ReadSNP} \link{ReturnUniqueBases}
#' @examples
#' data(fakeData)
#' splitfakeData <- SplitSNP(fakeData)
#' GetBaseFrequencies(splitfakeData[,1])
#' GetBaseFrequencies(splitfakeData[,2])
#' GetBaseFrequencies(splitfakeData[,3])

GetBaseFrequencies <- function(site){
  bases <- ReturnUniqueBases(site)
  return(sapply(bases, function(x) length(which(site == x))/length(site)))
}