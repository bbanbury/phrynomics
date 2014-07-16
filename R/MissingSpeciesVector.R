MissingSpeciesVector <- function(data, SpeciesNames){
#returns T/F vector, TRUE is good to use
  species <- GetSpecies(SpeciesNames)
  for(i in sequence(length(species))) {
    combinedLocus <- paste(data[grep(species[i], names(data))], collapse="")
    if(all(strsplit(combinedLocus, "")[[1]] == "N")) {
      return(FALSE)
    }
  }
  return(TRUE)
}
