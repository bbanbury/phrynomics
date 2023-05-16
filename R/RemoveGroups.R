#' Remove Groups
#' 
#' This function will remove taxon names that do not match with an ingroup ID.  For example, if all the ingroup taxa share a taxon flag identifier "PH", then this function will remove any that do not have this flag. Useful for comparing ingroup to whole data or for analyses like SNAPP that do not require ingroups.
#' @param SNPdataset SNP data in the class "matrix", "data.frame", or "snp"
#' @param groupFlag Unique flag for ingroup taxa
#' @export
#' @return Returns a subset dataset with only ingroup taxa
#' @seealso \link{ReadSNP} \link{WriteSNP}
#' @examples
#' data(fakeData)
#' RemoveGroups(fakeData, groupFlag ="in")

RemoveGroups <- function(SNPdataset, groupFlag=""){
  snpclass <- "table"
  if(inherits(SNPdataset, "snp")){
    snpclass <- "snp"
    SNPdataset <- SNPdataset$data
  }
  whichIngroup <- grep(groupFlag, rownames(SNPdataset))
  newData <- data.frame(SNPdataset[whichIngroup,], stringsAsFactors=FALSE)
  rownames(newData) <- rownames(SNPdataset)[whichIngroup]
  colnames(newData) <- colnames(SNPdataset)
  if(snpclass == "snp")
    return(ReadSNP(newData))
  else
    return(newData)
}