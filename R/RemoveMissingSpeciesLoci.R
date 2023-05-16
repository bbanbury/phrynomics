#' Remove Loci That Do Not have Species Coverage
#' 
#' This function will remove loci where entire species are missing data (across all sites of a locus). Species names are important for this function, so if they are read in incorrectly, if could affect the results. They should be in the format where each species shares a unique flag and are then numbered (for example, species1, species2, species3 would be three individuals of the same species). If you want to check species see the function \code{GetSpecies}. 
#' @param SNPdataset SNP data in the class "matrix", "data.frame", or "snp"
#' @param chatty Optional print to screen messages
#' @export
#' @return Returns a subset dataset with sites where at least one individual has data. 
#' @seealso \link{ReadSNP} \link{WriteSNP} \link{GetSpecies} \link{IsMissingSpecies}
#' @examples
#' data(fakeData)
#' RemoveMissingSpeciesLoci(fakeData)
#' RemoveMissingSpeciesLoci(fakeData, chatty=TRUE)

RemoveMissingSpeciesLoci <- function(SNPdataset, chatty=FALSE){
  snpclass <- "table"
  if (inherits(SNPdataset, "snp")){
    snpclass <- "snp"
    SNPdataset <- SNPdataset$data
  }
  speciesWithN <- apply(SNPdataset, 2, IsMissingSpecies, SpeciesNames = rownames(SNPdataset), chatty=chatty)
  if(chatty) 
    message(paste("removed", length(which(speciesWithN == "FALSE")), "of", length(speciesWithN), "SNPs"))
  newSNPdataset <- as.data.frame(SNPdataset[, speciesWithN])
  rownames(newSNPdataset) <- rownames(SNPdataset)
  if (snpclass == "snp") 
    return(ReadSNP(newSNPdataset))
  else return(newSNPdataset)
}
