#' Calculate Amount of Data Overlap
#' 
#' This function will calculate the amount of data overlap in a dataset.  This assesses how much data is present in a dataset and how that data might be distributed across taxa. 
#' @param SNPdataset SNP dataset in the class "matrix", "data.frame", or "snp"
#' @export
#' @return Returns a list of two items. The first item is a matrix of frequencies. The individual values are based on present/absent status. If a informational base (non-missing data) exists, then the frequency of present data is entered here.  For example, if all of the species in a site have data, then 1 will be entered for each one. Which base is present does not matter, only that there is or is not data present. This value is the same as the probability that the species contains data for that site.  The second item are the row/species means.  This value equates with the probability of an site overlapping with another species in the dataset. 
#' @seealso \link{ReadSNP} \link{MakePresentAbsent}
#' @examples
#' data(fakeData)
#' DataOverlap(fakeData)

DataOverlap <- function (SNPdataset){
  if(inherits(SNPdataset, "snp"))
    SNPdataset <- SNPdataset$data
  presAbsentDataset <- MakePresentAbsent(SNPdataset, calc="loci")
  rownames(presAbsentDataset) <- rownames(SNPdataset)
  totals <- apply(presAbsentDataset, 2, sum)
  probs <- totals/dim(SNPdataset)[1]
  for(i in sequence(dim(presAbsentDataset)[2])){  #probably could clean this up by using replace()
    presAbsentDataset[which(presAbsentDataset[,i] == 1),i] <- probs[i]
  }  
  speciesTotals <- apply(presAbsentDataset, 1, mean)
  return(list(presAbsentDataset=presAbsentDataset, speciesTotals=speciesTotals))
}
