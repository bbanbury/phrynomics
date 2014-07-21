#' Return Unique Bases
#' 
#' This function will return the unique nucleotide bases of an SNP accounting for an ambiguity code
#' @param SNP A single SNP
#' @export
#' @return Returns a character vector with base possibilities.
#' @seealso \link{ReadSNP} \link{WriteSNP} \link{ReturnNucs}
#' @examples
#' ReturnUniqueBases(c("A", "A", "R"))
#' 
#' data(fakeData)
#' fakeData <- SplitSNP(fakeData)
#' ReturnUniqueBases(fakeData[,1])
#' ReturnUniqueBases(fakeData[,12])

ReturnUniqueBases <- function(SNP){
  bases <- unique(c(sapply(SNP, ReturnNucs, forSNAPP=TRUE), recursive=T)) #forSNAPP arg
  if(any(bases == "-"))
    bases <- bases[-which(bases == "-")]
  return(bases)
}
