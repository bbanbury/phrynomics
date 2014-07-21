##' Make Data Present/Absent
#' 
#' This function will remove breaks between loci and change all data to boolean values.  If data is present in any base pair or ambiguous character, then the value is assigned present (1), if data is missing (ie, "N", "-", or "?") then it is assigned absent (0).   
#' @param SNPdataset SNP dataset in the class "matrix", "data.frame", or "snp"
#' @param missingChar The character that represents missing data, This function will search for missing data ("N", "-" or "?"). You can choose to leave the default "any" which will search for all three variables or select whichever are present in your dataset. 
#' @export
#' @return Returns a matrix with each column representing a single SNP with present/absent data. 
#' @seealso \link{ReadSNP} \link{ConvertMissingData}
#' @examples
#' data(fakeData)
#' MakePresentAbsent(fakeData)

MakePresentAbsent <- function(SNPdataset){
  if(class(SNPdataset) == "snp")
    SNPdataset <- SNPdataset$data
  split <- SplitSNP(SNPdataset)
  split <- split[,-which(split[1,] == " ")]
  split[c(which(split=="-"), which(split=="N"), which(split=="?"))] <- 0
  split[-c(which(split=="-"), which(split=="N"), which(split=="?"), which(split=="0"))] <- 1
  presabsentdata <- matrix(as.numeric(unlist(split)),nrow=nrow(split))
  rownames(presabsentdata) <- rownames(SNPdataset)
  return(presabsentdata)
}
