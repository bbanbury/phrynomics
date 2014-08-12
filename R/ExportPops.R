##' Export Population Allele Counts
#' 
#' This function will ...  
#' @param SNPdataset SNP data in the class "matrix", "data.frame", or "snp"
#' @param file A file name specified by either a variable or a double-quoted string. If not specified (ie, file=""), then  data is printed to the console. 
#' @param subsets A vector of flags which will subset the populations. 
#' @param includeFull logical, to include the full population data in addition to the subsets.  
#' @param ... additional arguments passed to \code{write.table}
#' @export
#' @return Returns a table with allele counts for a whole population and any subset populations. 
#' @seealso \link{ReadSNP} \link{AddAFlag} \link{ReturnAlleleCounts}
#' @examples
#' data(fakeData)
#' #fakeData <- RemoveNonBinary(RemoveInvariantSites(fakeData))
#' #fakeData <- AddAFlag(fakeData, flagToAdd="WEST", taxa=c("in1", "in2", "in3"))
#' #ExportPops(fakeData, subsets=c("WEST"))

ExportPops <- function(SNPdataset, file="", subsets=c(), includeFull=FALSE, ...){
#this will print to screen and write to file a treemix file.
#subsets should be the flags ex: c("north", "west").  Will write full pop then n subsets
  if(class(SNPdataset) == "snp")
    SNPdataset <- SNPdataset$data
  alleleCounts <- lapply(SplitSNP(cSNP(SNPdataset)), ReturnAlleleCounts)
  pops <- matrix(ncol=1+length(subsets), nrow=length(alleleCounts))
  newRownames <- NULL
  for(i in sequence(dim(pops)[1])){
    pops[i,1] <- paste(alleleCounts[[i]], collapse=",")
    newRownames <- c(newRownames, paste(names(alleleCounts[[i]]), collapse=","))
  }
  rownames(pops) <- newRownames
  colnames(pops) <- c("full", subsets)
    for(j in sequence(length(subsets))){
      subsetPop <- RemoveGroups(SNPdataset, groupFlag=subsets[j])
      alleleCounts <- lapply(SplitSNP(cSNP(subsetPop)), ReturnAlleleCounts)
      alleleCounts <- addazero(alleleCounts, rownames(pops))
      for(row in sequence(dim(pops)[1])){
        pops[row,j+1] <- paste(alleleCounts[[row]], collapse=",")
      }
    }
  if(!includeFull)
    pops <- pops[,-1]  
  if(!is.null(file))
    write.table(pops, file=file, row.names=FALSE, quote=FALSE, ...)
  return(pops)
}
