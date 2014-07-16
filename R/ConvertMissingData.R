##' Convert Missing Data Character
#' 
#' This function will convert any character in a dataset to another.  This is helpful when one program outputs "N" representing missing data and the next requires "-".  You could also use this for converting other types of characters too (for example, if you wanted "A" to be "O"). 
#' @param SNPdataset SNP data in the class "matrix", "data.frame", or "snp"
#' @param oldmissing The character that represents missing data
#' @param newmissing The new character to represent missing data
#' @export
#' @return Returns a dataset of the same dimentions as the original, but with swapped out characters. 
#' @seealso \link{ReadSNP} \link{CalculateMissingData}
#' @examples
#' data(fakeData)
#' ConvertMissingData(fakeData)

ConvertMissingData <- function(SNPdataset, oldmissing="N", newmissing="?"){
  if(class(SNPdataset) == "snp")
    SNPdataset <- SNPdataset$data
  for(i in sequence(dim(SNPdataset)[1])){
    SNPdataset[i,] <- gsub(oldmissing, newmissing, SNPdataset[i,])
  }
  return(SNPdataset)
}
