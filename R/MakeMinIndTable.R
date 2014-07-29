MakeMinIndTable <- function(SNPdataset, showEvery=5, missingChar="any"){
  if(class(SNPdataset) == "snp")
    SNPdataset <- SNPdataset$data
  toKeep <- rev(c(1, seq.int(to=dim(SNPdataset)[1], from=showEvery, by=showEvery)))
  indTable <- matrix(nrow=length(toKeep), ncol=4)
  rownames(indTable) <- paste("c", toKeep, sep="")
  colnames(indTable) <- c("MinIndivduals", "Loci", "VariableSites", "MissingLoci(%)")
  indTable[,1] <- toKeep
  sumloci <- apply(MakePresentAbsent(SNPdataset, calc="loci", missingChar=missingChar), 2, sum)
  sumVarsites <- apply(MakePresentAbsent(RemoveInvariantSites(SNPdataset), calc="sites", missingChar=missingChar), 2, sum)
  for(row.ind in sequence(dim(indTable)[1])){
    a <- length(which(sumloci >= indTable[row.ind, 1]))
    b <- length(which(sumVarsites >= indTable[row.ind, 1]))
    c <- mean(CalculateMissingData(ReduceMinInd(SNPdataset, calc="loci", indTable[row.ind,1]), "loci"))
    indTable[row.ind, c(2,3,4)] <- c(a,b,c)
  }
  return(indTable)
}
#MakeMinIndTable(fakeData, showEvery=1)
#MakeMinIndTable(l, showEvery=5)