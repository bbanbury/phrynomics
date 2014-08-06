#' Add A Zero to alleleCounts list
#' 
#' This function will add a zero if allele counts are missing from subpops
#' @param SNPdataset SNP data in the class "matrix", "data.frame", or "snp"
#' @param rownamesOriginal rownames from the original full dataset.  These will be used to compare the subset to. 
#' @export
#' @return Returns a new alleleCounts list with zeros in the right position
#' @seealso \link{ReadSNP} \link{WriteSNP} \link{ExportPops}

addazero <- function(alleleCounts, rownamesOriginal){
#this function will add missing alleles to subset data when the smaller population doesn't have the allele present.  For example, if the total pop has A/G, but the subset is only A, then the allele counts will only return 100A rather than 100A/0G.  
#alleleCounts should be from the function ReturnAlleleCounts 
#rownamesOriginal should be from pops, but would also be "A,G" for each locus
  for(i in sequence(length(alleleCounts))){
    if(paste(names(alleleCounts[[i]]), collapse=",") != rownamesOriginal[i]){
      if(length(alleleCounts[[i]]) == 0){  #no data
        alleleCounts[[i]] <- c(0,0)
        names(alleleCounts[[i]]) <- strsplit(rownamesOriginal[i], ",")[[1]]
      }
      if(length(grep("[,]", paste(names(alleleCounts[[i]]), collapse=","))) == 0){  #no comma means one allele
        position <- grep(paste(names(alleleCounts[[i]]), collapse=","), strsplit(rownamesOriginal[i], ",")[[1]])
        if(position == 1){
          alleleCounts[[i]] <- c(alleleCounts[[i]], 0)
          names(alleleCounts[[i]])[2] <- strsplit(rownamesOriginal[i], ",")[[1]][2]
        }
        if(position == 2){
          alleleCounts[[i]] <- c(0, alleleCounts[[i]])
          names(alleleCounts[[i]])[1] <- strsplit(rownamesOriginal[i], ",")[[1]][1]
        }
      }
      if(length(grep("[,]", paste(names(alleleCounts[[i]]), collapse=","))) == 1){  #one comma means two alleles
        alleleCounts[[i]] <- rev(alleleCounts[[i]])
      }

    }
  }
  return(alleleCounts)
}
