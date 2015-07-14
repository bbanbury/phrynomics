#' Missing Species Vector
#' 
#' This function will determine if entire species are missing data (across all sites of a locus). Species names are important for this function, so if they are read in incorrectly, if could affect the results. They should be in the format where each species shares a unique flag and are then numbered (for example, species1, species2, species3 would be three individuals of the same species). If you want to check species see the function \code{GetSpecies}. 
#' @param locus A single locus (can have multiple sites)
#' @param SpeciesNames Vector of species names that will cluster individuals. This will likely be rownames(SNPdataset)
#' @param chatty Option to print details to screen
#' @export
#' @return Returns a TRUE/FALSE vector for each locus in the dataset 
#' @seealso \link{ReadSNP} \link{WriteSNP} \link{GetSpecies} \link{RemoveMissingSpeciesLoci}
#' @examples
#' data(fakeData)
#' Spnames <- rownames(fakeData)
#' IsMissingSpecies(fakeData[,1], Spnames)
#' IsMissingSpecies(fakeData[,2], Spnames)
#' IsMissingSpecies(fakeData[,3], Spnames)

IsMissingSpecies <- function(locus, SpeciesNames, chatty=FALSE){
  species <- GetSpecies(SpeciesNames)
  for(i in sequence(length(species))) {
    combinedLocus <- paste(locus[grep(species[i], SpeciesNames)], collapse="")
    if(all(strsplit(combinedLocus, "")[[1]] == "N")) {
      if(chatty)
        print(paste("Species", species[i], "has all missing"))
      return(FALSE)
    }
  }
  return(TRUE)
}
