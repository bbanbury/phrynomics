#' Remove Invariant Sites
#' 
#' This function will determine if a site is variable (TRUE) or not (FALSE) and remove those sites that are not variable. Ambiguity codes are taken into account to mean either heterozygotes or uncertainty. For example, c("A", "A", "S") will return TRUE, because "S" reduces to "G" and "C" and so either of those is variable with "A"; c("A", "A", "M") will return FALSE, because M can be either A or C. If the "M" is uncertain and is an "A" then it is not variable.  
#' @param SNPdataset SNP data in the class "matrix", "data.frame", or "snp"
#' @param chatty Optional print to screen messages
#' @export
#' @return Returns a subset dataset with only variable sites (ie SNPs). 
#' @seealso \link{ReadSNP} \link{WriteSNP} \link{IsVariable}
#' @examples
#' data(fakeData)
#' RemoveInvariantSites(fakeData)
#' RemoveInvariantSites(fakeData, chatty=TRUE)

RemoveInvariantSites <- function(SNPdataset, chatty=FALSE){
  snpclass <- "table"
  if(class(SNPdataset) == "snp"){
    snpclass <- "snp"
    SNPdataset <- SNPdataset$data
  }
  snps <- sum(nchar(SNPdataset[1,]))
  initialLociLengths <- nchar(SNPdataset[1,])
  splitdata <- SplitSNP(SNPdataset)
  KeepVector <- apply(splitdata, 2, IsVariable)  
  breaks <- which(splitdata[1,] == " ")
  newSNPdataset <- cSNP(splitdata, KeepVector=KeepVector, maintainLoci=TRUE)
  newsnps <- sum(nchar(newSNPdataset[1,]))
  if(chatty)
    message(paste("removed", snps-newsnps, "of", snps, "sites"))
  if(snpclass == "snp")
    return(ReadSNP(newSNPdataset))
  else
    return(newSNPdataset)
}
