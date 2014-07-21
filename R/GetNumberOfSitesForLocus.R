##' Get The Number SNPs in a Locus
#' 
#' This function will calculate The number of SNPs for a given locus. 
#' @param SNPdataset SNP dataset in the class "matrix", "data.frame", or "snp"
#' @param locus Which locus to calculate
#' @export
#' @return Returns a number of SNPs 
#' @seealso \link{ReadSNP} \link{TakeSingleSNPfromEachLocus}
#' @examples
#' data(fakeData)
#' GetNumberOfSitesForLocus(fakeData, 1)
#' GetNumberOfSitesForLocus(fakeData, 2)
#' GetNumberOfSitesForLocus(fakeData, 3)

GetNumberOfSitesForLocus <- function(SNPdataset, locus){
  if(class(SNPdataset) == "snp")
    SNPdataset <- SNPdataset$data
  return(nchar(SNPdataset[1,])[locus])
}
