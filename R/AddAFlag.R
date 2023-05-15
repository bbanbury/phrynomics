#' Add A Flag To A Subset Group
#' 
#' This function will add a flag to the names taxa in the list. 
#' @param SNPdataset SNP data in the class "matrix", "data.frame", or "snp"
#' @param flagToAdd Unique flag for ingroup taxa
#' @param taxa Vector of names
#' @export
#' @return Returns a subset dataset with only ingroup taxa
#' @seealso \link{ReadSNP} \link{WriteSNP} \link{ExportPops}
#' @examples
#' data(fakeData)
#' AddAFlag(fakeData, flagToAdd="WEST", taxa=c("in1", "in2", "in3"))

AddAFlag <- function(SNPdataset, flagToAdd="", taxa=c()){
#will add a flag to certain taxa for subsetting later
  snpclass <- "table"
  if(inherits(SNPdataset, "snp")){
    snpclass <- "snp"
    SNPdataset <- SNPdataset$data
  }
  whichTaxa <- which(rownames(SNPdataset) %in% taxa)
  rownames(SNPdataset)[whichTaxa] <- paste(flagToAdd, rownames(SNPdataset)[whichTaxa], sep="")
  if(snpclass == "snp")
    return(ReadSNP(SNPdataset))
  else
    return(SNPdataset)
}
