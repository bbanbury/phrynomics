##' Return Allele COunts
#' 
#' This function will calculate the number of alleles in a population. Heterozygotes get split between their bases.   
#' @param site Single SNP site as a vector
#' @export
#' @return Returns a named list of allele counts.
#' @seealso \link{ReadSNP} \link{GetBaseFrequencies} \link{ExportPops}
#' @examples
#' data(fakeData)
#' splitfakeData <- SplitSNP(fakeData)
#' ReturnAlleleCounts(splitfakeData[,1])
#' ReturnAlleleCounts(splitfakeData[,2])
#' lapply(splitfakeData, ReturnAlleleCounts)

ReturnAlleleCounts <- function(site){
  uniques <- unique(site)
  bases <- ReturnUniqueBases(site)
  allelesToCount <- NULL
  for(i in sequence(length(uniques))){
    if(all(ReturnNucs(uniques[i]) %in% bases))
      allelesToCount <- c(allelesToCount, uniques[i])
  }
  counts <- sapply(allelesToCount, function(x) length(which(site == x)))
  whichToDouble <- which(names(counts) %in% c("A", "C", "T", "G"))
  for(i in whichToDouble){   #double homos
    counts[i] <- counts[i] * 2
  }
  for(i in sequence(length(counts))){
    if(!(i %in% whichToDouble)){
      counts <- counts + counts[i]  #add a het to each other allele
      counts <- counts[-i]
    }
  }
  return(counts)
}
