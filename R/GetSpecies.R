GetSpecies <- function(taxa){
  return(unique(gsub("\\d$", "", taxa)))
}
