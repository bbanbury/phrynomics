#' Percent Missing Data
#' 
#' This function will determine the percent missing from a locus.    
#' @param locus A single locus
#' @param missingChar The character that represents missing data, This function will search for missing data ("N", "-" or "?"). You can choose to leave the default "any" which will search for all three variables or select whichever are present in your dataset. 
#' @export
#' @return numerical, A single value for percentage of missing data (0 is none, 1 is all) 
#' @seealso \link{ReadSNP} \link{MakePresentAbsent} \link{CalculateMissingData} \link{IsMissing}
#' @examples
#' data(fakeData)
#' PercentMissing(fakeData[1,1])
#' IsMissing(fakeData[1,3])

PercentMissing <- function(locus, missingChar="any"){
  if (missingChar == "any")
    missingChar <- c("N", "-", "?")
  char <- match.arg(arg=missingChar, choices=c("N", "-", "?"), several.ok=TRUE)  
  if(length(grep("[N?-]", locus)) > 0) {
    return(length(grep("[N?-]", strsplit(locus, "")[[1]]))/length(strsplit(locus, "")[[1]]))
  }
  else
    return(0)
}