#' Reduce Data by Minimum Individuals
#' 
#' This function will reduce a dataset to include sites or loci that meet the minimum number of individuals.  For example, if you want only data for which at least 75% of the individuals in your dataset are represented, then it will remove those loci/sites that have less.     
#' @param SNPdataset SNP dataset in the class "matrix", "data.frame", or "snp"
#' @param showEvery An integer to describe how many datasets to calculate.  If \code{showEvery=3}, you will get a table that returns minimum individual data increasing by 3 individuals each time (ex: 1, 4, 7, 10, etc...)
#' @param missingChar The character that represents missing data, This function will search for missing data ("N", "-" or "?"). You can choose to leave the default "any" which will search for all three variables or select whichever are present in your dataset. 
#' @param file An optional agument to write the table to a file. 
#' @param ... Optional additional arguments to pass to file writing (such as appending the file with more data).   
#' @export
#' @return Returns a matrix with the minimum number of individuals, the number of loci, variable sites, and percent of missing data for that number of individuals.  
#' @seealso \link{ReadSNP} \link{ReduceMinInd} \link{plotMissing}
#' @examples
#' data(fakeData)
#' MakeMinIndTable(fakeData, showEvery=1)
#' MakeMinIndTable(fakeData, showEvery=3)

MakeMinIndTable <- function(SNPdataset, showEvery=5, missingChar="any", file=NULL, ...){
  if(inherits(SNPdataset, "snp"))
    SNPdataset <- SNPdataset$data
  toKeep <- rev(c(1, seq.int(to=dim(SNPdataset)[1], from=showEvery, by=showEvery)))
  indTable <- matrix(nrow=length(toKeep), ncol=4)
  rownames(indTable) <- paste("c", toKeep, sep="")
  colnames(indTable) <- c("MinIndivduals", "Loci", "VariableSites", "MissingLoci(%)")
  indTable[,1] <- toKeep
  sumloci <- apply(MakePresentAbsent(SNPdataset, calc="loci", missingChar=missingChar), 2, sum)
  sumVarsites <- apply(MakePresentAbsent(RemoveInvariantSites(SNPdataset), calc="sites", missingChar=missingChar), 2, sum)
  for(row.ind in sequence(dim(indTable)[1])){
    SNPdataset <- ReduceMinInd(SNPdataset, calc="loci", indTable[row.ind,1])
    a <- length(which(sumloci >= indTable[row.ind, 1]))
    b <- length(which(sumVarsites >= indTable[row.ind, 1]))
    c <- mean(CalculateMissingData(SNPdataset, "loci"))
    indTable[row.ind, c(2,3,4)] <- c(a,b,c)
  }
  if(!is.null(file))
    write.table(indTable, file=file, ...)
  return(indTable)
}