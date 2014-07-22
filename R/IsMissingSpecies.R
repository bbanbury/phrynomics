#' Missing Species Vector
#' 
#' This function will determine if entire species are missing data (across all SNPs of a locus). Species names are important for this function, so if they are read in incorrectly, if could affect the results. They should be in the format where each species shares a unique flag and are then numbered (for example, species1, species2, species3 would be three individuals of the same species). If you want to check species see the function \code{GetSpecies}. 
#' @param SNPlocus A single SNP locus (can have multiple SNPs)
#' @param SpeciesNames Vector of species names that will cluster individuals. This will likely be rownames(SNPdataset)
#' @export
#' @return Returns a TRUE/FALSE vector for each locus in the dataset 
#' @seealso \link{ReadSNP} \link{WriteSNP} \link{GetSpecies} \link{RemoveMissingSpeciesLoci}
#' @examples
#' data(fakeData)
#' Spnames <- rownames(fakeData)
#' IsMissingSpecies(fakeData[,1], Spnames)
#' IsMissingSpecies(fakeData[,2], Spnames)
#' IsMissingSpecies(fakeData[,3], Spnames)

IsMissingSpecies <- function(SNPlocus, SpeciesNames){
  species <- GetSpecies(SpeciesNames)
  for(i in sequence(length(species))) {
    combinedLocus <- paste(SNPlocus[grep(species[i], names(SNPlocus))], collapse="")
    if(all(strsplit(combinedLocus, "")[[1]] == "N")) {
      return(FALSE)
    }
  }
  return(TRUE)
}
#fix plotting function
#add third panel with missing data info
#add binary/variable site data for plot
#change "snp" to "sites"
#proofread in textwrangler
#create package
#test/check package
#check shiny site using new functions
