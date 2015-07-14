##' Get Species Names in a Dataset
#' 
#' This function will detect the names of species in a dataset by removing the number identifyer.  
#' @param taxa A character vector of taxon names. Numbers in species names represent individuals of the same species, for example "taxon1" and "taxon2" are both the same species "taxon".  
#' @export
#' @return Returns a vector of unique species names
#' @seealso \link{ReadSNP} \link{IsMissingSpecies} \link{RemoveMissingSpeciesLoci}
#' @examples
#' taxa <- c("taxon1", "taxon2", "taxon3")
#' GetSpecies(taxa)
#' 
#' data(fakeData)
#' GetSpecies(rownames(fakeData))

GetSpecies <- function(taxa){
  return(unique(gsub("\\d*$", "", taxa)))
}
