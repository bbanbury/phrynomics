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
#' ConvertMissingData(fakeData, "N", "O")
#' ConvertMissingData(fakeData, c("N", "-"), "?")
#' ConvertMissingData(fakeData, "any", "?")

ConvertMissingData <- function(SNPdataset, oldmissing=c("-","N"), newmissing="N"){
  snpclass <- "table"
  if(inherits(SNPdataset, "snp")){
    snpclass <- "snp"
    SNPdataset <- SNPdataset$data
  }
  if(any(oldmissing == "any"))
    oldmissing <- c("N", "-", "\\?")
  char <- match.arg(arg=oldmissing, choices=c("N", "-", "\\?"), several.ok=TRUE)
  pattern <- paste("[", paste(char, collapse=""), "]", sep="")
  newSNPdataset <- data.frame(lapply(SplitSNP(SNPdataset), gsub, pattern=pattern, replacement=newmissing), stringsAsFactors=FALSE)
  newSNPdataset <- cSNP(newSNPdataset)
  rownames(newSNPdataset) <- rownames(SNPdataset)
  colnames(newSNPdataset) <- colnames(SNPdataset)
  if(snpclass == "snp")
    return(ReadSNP(newSNPdataset))
  else
    return(newSNPdataset)
}
