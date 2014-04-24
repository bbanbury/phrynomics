TranslateBases <- function(data, translateMissing=TRUE, translateMissingChar="-", ordered=FALSE){
#here "ordered" characters means that per SNP, 0 will be the most frequent base, 1 will be the heterozygote, and 2 will be the less frequent base.  Used for SNAPP.  
#If ordered=FALSE, then A will always return a 1, T=2, G=3, C=4, and amby codes get parentheses.  Used for MrBayes.  
  catData <- cSNP(data, maintainLoci=FALSE)
  splitdata <- SplitSNP(catData)
  if(translateMissing){
    splitdata[which(splitdata == "N")] <- translateMissingChar  #translate missing data
    splitdata[which(splitdata == "-")] <- translateMissingChar  
    splitdata[which(splitdata == "?")] <- translateMissingChar
  }
  if(ordered){
    bases <- apply(splitdata, 2, ReturnUniqueBases)
    baseFrequencies <- apply(splitdata, 2, GetBaseFrequencies)
    for(i in sequence(dim(splitdata)[2])){
      if(baseFrequencies[1,i] == baseFrequencies[2,i]){
        MostFrequentBase <- bases[1,i]
        LeastFrequentBase <- bases[2,i]       
      }
      if(baseFrequencies[1,i] != baseFrequencies[2,i]){
        MostFrequentBase <- bases[which(baseFrequencies[,i] == max(baseFrequencies[,i])), i]
        LeastFrequentBase <- bases[which(baseFrequencies[,i] == min(baseFrequencies[,i])), i]
      }  
      zeros <- which(splitdata[,i] == MostFrequentBase)
      splitdata[zeros,i] <- 0
      if(any(splitdata[,i] == ReturnAmbyCode(bases[,i]))){
        ones <- which(splitdata[,i] == ReturnAmbyCode(bases[,i]))
        splitdata[ones,i] <- 1
      }
      if(any(which(splitdata[,i] == LeastFrequentBase))) {
        twos <- which(splitdata[,i] == LeastFrequentBase)
        splitdata[twos,i] <- 2
      }
    }    
  }
  else {
    for(i in sequence(dim(splitdata)[1])){
      splitdata[i,] <- sapply(splitdata[i,], ReturnMrBayesAmbyCode)
    }
  }
  return(cSNP(splitdata))
}
