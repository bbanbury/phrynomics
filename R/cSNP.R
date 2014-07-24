#' Combine SNP
#' 
#' Combine split character strings into single character elements. For example, the string "A" "T" "A" will return "ATA". Breaks in loci are maintained using \code{maintainLoci=TRUE}.
#' @param splitSNP SNP data matrix in the class "data.frame", "matrix", or "snp", where each site is a separate column
#' @param KeepVector A boolean vector for returning specific sites.  If left to NULL, all SNPs will be returned. 
#' @param maintainLoci If \code{maintainLoci=TRUE}, loci will be split by a space in the returned matrix. 
#' @export
#' @return Returns a data frame of sites either with or without spaces indicating breaks in loci.  
#' @seealso \link{ReadSNP} \link{SplitSNP}
#' @examples
#' data(fakeData)
#' fakeData <- ReadSNP(fakeData)
#' cSNP(SplitSNP(fakeData))
#' cSNP(SplitSNP(fakeData), maintainLoci=FALSE)
#' 
#' coolData <- matrix(data=c("A", "T", "C", "A", "S", "N"), nrow=2)
#' cSNP(coolData)

cSNP <- function(splitSNP, KeepVector=NULL, maintainLoci=TRUE){
  snpclass <- "table"
  if(class(splitSNP) == "snp"){
      snpclass <- "snp"
      splitSNP <- splitSNP$data
  }
   if(is.null(KeepVector))
      KeepVector <- rep(TRUE, dim(splitSNP)[2])
  catSNP <- matrix(nrow=dim(splitSNP)[1], ncol=length(which(splitSNP[1,] == " "))+1) 
  rownames(catSNP) <- rownames(splitSNP)
  for(j in sequence(dim(catSNP)[1])) {
    #totally annoying strsplit issue where it cleaves off one space at the end
    string <- strsplit(paste(splitSNP[j,][which(KeepVector=="TRUE")], collapse=""), " ")[[1]]
    if(length(string) +1 == dim(catSNP)[2]){
      string <- c(string, "")
    }
    catSNP[j,] <- string
  }
  if(any(which(catSNP[1,] == "")))
    catSNP <- catSNP[,-which(catSNP[1,] == "")]  #rewrite with lost loci
  if(!maintainLoci)
    catSNP <- apply(catSNP, 1, paste, collapse="")
  if(snpclass == "snp")
    return(ReadSNP(data.frame(catSNP, stringsAsFactors=FALSE)))
  else
    return(data.frame(catSNP, stringsAsFactors=FALSE))
}