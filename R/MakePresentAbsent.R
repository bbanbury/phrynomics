#' Make Data Present/Absent
#' 
#' This function will remove breaks between loci and change all data to boolean values.  If data is present in any base pair or ambiguous character, then the value is assigned present (1), if data is missing (ie, "N", "-", or "?") then it is assigned absent (0).   
#' @param SNPdataset SNP dataset in the class "matrix", "data.frame", or "snp"
#' @param calc Calculate missing data either by individual "sites" or "loci". If you choose to calculate missing data by loci, any site with information will count (for example, "NNNNA" will count as non-missing)
#' @param missingChar The character that represents missing data, This function will search for missing data ("N", "-" or "?"). You can choose to leave the default "any" which will search for all three variables or select whichever are present in your dataset. 
#' @export
#' @return Returns a matrix with each column representing a single site with present/absent data. 
#' @seealso \link{ReadSNP} \link{ConvertMissingData} \link{CalculateMissingData}
#' @examples
#' data(fakeData)
#' MakePresentAbsent(fakeData, "sites")
#' MakePresentAbsent(fakeData, "loci")

MakePresentAbsent <- function(SNPdataset, calc="loci", missingChar="any"){
  if(inherits(SNPdataset, "snp"))
    SNPdataset <- SNPdataset$data
  if(calc == "sites"){
    SNPdataset <- SplitSNP(SNPdataset)
    if(length(which(SNPdataset[1,] == " ")) > 0)
      SNPdataset <- SNPdataset[,-which(SNPdataset[1,] == " ")]
  }
  presabsentdata <- apply(SNPdataset, c(1,2), IsMissing, missingChar=missingChar)
  return(presabsentdata)
}