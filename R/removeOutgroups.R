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
#' removeGroup(fakeData, taxon="in")

removeGroup <- function(SNPdataset, groupFlag=""){
  if(class(SNPdataset) == "snp"){
    newData <- SNPdataset$data[grep(groupFlag, rownames(SNPdataset$data)),]
    newData <- ReadSNP(newData)
    return(newData)
  }
  else
    return(SNPdataset[grep(groupFlag, rownames(SNPdataset)),])
}