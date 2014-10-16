##' Return Allele COunts
#' 
#' This function will calculate the number of alleles in a population. Heterozygotes get split between their bases.   
#' @param site Single SNP site as a vector
#' @param allelenames whether to return a named vector
#' @param bgc For Bayesian Genomic Clines (Gompert and Buerkle 2011)
#' @export
#' @return Returns a named list of allele counts.
#' @seealso \link{ReadSNP} \link{GetBaseFrequencies} \link{ExportPops}
#' @examples
#' data(fakeData)
#' splitfakeData <- SplitSNP(fakeData)
#' ReturnAlleleCounts(splitfakeData[,1])
#' ReturnAlleleCounts(splitfakeData[,2])
#' lapply(splitfakeData, ReturnAlleleCounts)

ReturnAlleleCounts <- function(site, allelenames=TRUE, bgc=FALSE){
  uniques <- unique(site)
  bases <- ReturnUniqueBases(site)
  if(bgc){
    if(length(bases) == 0)
      return(c(-9, -9))
  }
  if(!bgc){
    if(length(bases) == 0)
      return(c(0, 0))
  }
  counts <- c(0,0)
  for(i in sequence(length(bases))){
    counts[i] <- counts[i] + length(which(site == bases[i])) *2
  }
  if(any(!uniques %in% bases)){
    whichambys <- ReturnAmbyCode(bases)
    counts <- counts + length(which(site == uniques[whichambys]))
  }
  names(counts) <- bases
  if(!bgc){
    if(any(is.na(names(counts))))
      counts <- counts[-which(is.na(names(counts)))]
  }
  return(counts)
}