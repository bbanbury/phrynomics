#' Remove Non-binary Sites
#' 
#' This function will determine if a site is binary (TRUE) or not (FALSE) and remove those sites that are not binary.  For example, c("A", "A", "T") will be kept; c("A", "A", "A") or c("A", "G", "T") will be removed. Ambiguity codes are taken into account to mean either heterozygotes or uncertainty. 
#' @param SNPdataset SNP data in the class "matrix", "data.frame", or "snp"
#' @param chatty Optional print to screen messages
#' @export
#' @return Returns a subset dataset with only variable sites. 
#' @seealso \link{ReadSNP} \link{WriteSNP} \link{RemoveInvariantSites} \link{IsBinary}
#' @examples
#' data(fakeData)
#' RemoveNonBinary(fakeData)
#' RemoveNonBinary(fakeData, chatty=TRUE)

RemoveNonBinary <- function(SNPdataset, chatty=FALSE){
  snpclass <- "table"
  if(inherits(SNPdataset, "snp")){
    snpclass <- "snp"
    SNPdataset <- SNPdataset$data
  }
  snps <- sum(nchar(SNPdataset[1,]))
  splitdata <- SplitSNP(SNPdataset)
  binaryVector <- apply(splitdata, 2, IsBinary)
  newSNPdataset <- cSNP(splitdata, KeepVector=binaryVector)
  newsnps <- sum(nchar(newSNPdataset[1,]))
  if(chatty)
    message(paste("removed", snps-newsnps, "of", snps, "sites"))
  if(snpclass == "snp")
    return(ReadSNP(newSNPdataset))
  else
  return(newSNPdataset)
}
