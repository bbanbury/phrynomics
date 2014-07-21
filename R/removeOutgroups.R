#' Remove Outgroups
#' 
#' This function will remove taxon names that do not match with an ingroup ID.  FOr example, if all the ingroup taxa share a taxon flag identifier "PH", then this function will remove any that do not have this flag. Useful for comparing ingroup to whole tree data or for analyses like SNAPP that do not require ingroups.
#' @param SNPdataset SNP data in the class "matrix", "data.frame", or "snp"
#' @param taxon Unique flag for ingroup taxa
#' @export
#' @return Returns a subset dataset with only ingroup taxa
#' @seealso \link{ReadSNP} \link{WriteSNP}
#' @examples
#' data(fakeData)
#' removeOutgroups(fakeData, taxon="in")

removeOutgroups <- function(SNPdataset, taxon="PH"){
  if(class(SNPdataset) == "snp"){
    SNPdataset$data <- SNPdataset$data[grep(taxon, rownames(SNPdataset$data)),]
    return(SNPdataset)
  }
  return(SNPdataset[grep(taxon, rownames(SNPdataset)),])
}