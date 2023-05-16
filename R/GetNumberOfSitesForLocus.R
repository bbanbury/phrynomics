##' Get The Number Sites in a Locus
#' 
#' This function will calculate The number of sites for a given locus. 
#' @param SNPdataset SNP dataset in the class "matrix", "data.frame", or "snp"
#' @param locus Which locus to calculate
#' @export
#' @return Returns a number of sites 
#' @seealso \link{ReadSNP} \link{TakeSingleSNPfromEachLocus}
#' @examples
#' data(fakeData)
#' GetNumberOfSitesForLocus(fakeData, 1)
#' GetNumberOfSitesForLocus(fakeData, 2)
#' GetNumberOfSitesForLocus(fakeData, 3)

GetNumberOfSitesForLocus <- function(SNPdataset, locus){
  if(inherits(SNPdataset, "snp"))
    SNPdataset <- SNPdataset$data
  return(nchar(SNPdataset[1,])[locus])
}
