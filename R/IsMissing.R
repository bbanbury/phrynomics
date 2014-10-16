#' Is Missing Data
#' 
#' This function will determine if data is all missing from a locus.    
#' @param locus A single locus
#' @param missingChar The character that represents missing data, This function will search for missing data ("N", "-" or "?"). You can choose to leave the default "any" which will search for all three variables or select whichever are present in your dataset. 
#' @export
#' @return logical, return of whether is all missing data. For example, "AAAN" will be FALSE, "NNNN" will be TRUE. 
#' @seealso \link{ReadSNP} \link{MakePresentAbsent} \link{CalculateMissingData}
#' @examples
#' data(fakeData)
#' IsMissing(fakeData[1,1])
#' IsMissing(fakeData[1,3])

IsMissing <- function(locus, missingChar="any"){
  if (missingChar == "any")
    missingChar <- c("N", "-", "?")
  char <- match.arg(arg=missingChar, choices=c("N", "-", "?"), several.ok=TRUE)  
  if(length(grep("[N?-]", locus)) > 0 && all(unique(strsplit(locus, "")[[1]]) %in% char))
      return(0) 
  else
    return(1)
}