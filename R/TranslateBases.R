#' Translate Bases
#' 
#' This function will translate the bases from alpha characters to numeric. Several options exist for transformation, such as transforming missing characters and ordering characters by frequency. 
#' @param SNPdataset SNP data in the class "matrix", "data.frame", or "snp"
#' @param translateMissing Logical. Will transform missing characters if TRUE
#' @param translateMissingChar The character missing data will transform to. This has limited functionality, in that any missing data character (ie, "N", "-", and "?") will be transformed to the named character. For more options see \code{ConvertMissingData}.
#' @param ordered Logical. If ordered = FALSE (default), Sites will be transformed where A=1, T=2, G=3, and C=4. In the case of ambiguous characters, the possible bases are returned in parentheses. This is the format for MrBayes MkV models. If ordered = TRUE, each site is evaluated for base frequencies and the base with the most frequent occurence is assigned a 0 and the less frequent base is assigned a 1.  THis only works for binary traits.  Any heterozygotes (ambiguity with the present two bases) are assigned a 1.  This is the format for SNAPP analyses.  
#' @export
#' @return Returns a translated dataset.
#' @seealso \link{ReadSNP} \link{WriteSNP} \link{RemoveNonBinary}
#' @examples
#' data(fakeData)
#'
#' # Translate to MrBayes MkV formatting
#' TranslateBases(fakeData, translateMissing=FALSE, ordered=FALSE)
#' TranslateBases(fakeData, translateMissing=TRUE, translateMissingChar="?")
#' 
#' # Translate to SNAPP formatting
#' fakeData <- RemoveNonBinary(fakeData)
#' TranslateBases(fakeData, translateMissing=TRUE, ordered=TRUE)

TranslateBases <- function(SNPdataset, translateMissing=TRUE, translateMissingChar="-", ordered=FALSE){
  snpclass <- "table"
  if(class(SNPdataset) == "snp"){
    snpclass <- "snp"
    SNPdataset <- SNPdataset$data
  }
  catData <- cSNP(SNPdataset, maintainLoci=FALSE)
  splitdata <- as.matrix(SplitSNP(catData))
  if(translateMissing){
    splitdata <- ConvertMissingData(splitdata, "-", translateMissingChar)  #translate missing data
    splitdata <- ConvertMissingData(splitdata, "N", translateMissingChar)  #translate missing data
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
  if(snpclass == "snp")
    return(ReadSNP(cSNP(splitdata)))
  else
    return(cSNP(splitdata))
}
