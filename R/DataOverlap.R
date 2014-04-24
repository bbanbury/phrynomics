DataOverlap <- function (dataset){
#This function will calculate the amount of shared SNPs
  presAbsentDataset <- MakePresentAbsent(dataset)
  rownames(presAbsentDataset) <- rownames(dataset)
  totals <- apply(presAbsentDataset, 2, sum)
  probs <- totals/73
  for(i in sequence(dim(presAbsentDataset)[2])){  #probably could clean this up by using replace()
    presAbsentDataset[which(presAbsentDataset[,i] == 1),i] <- probs[i]
  }  
  speciesTotals <- apply(presAbsentDataset, 1, mean)
  return(list(presAbsentDataset, speciesTotals))
}
