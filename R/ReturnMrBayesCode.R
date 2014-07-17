#' Return MrBayes Codes
#' 
#' This function will take a nucleotide base and return a corresponding number
#' @param base A nucleotide base
#' @export
#' @return Returns a single number that corresponds with a base pair
#' @seealso \link{ReadSNP} \link{WriteSNP} \link{ReturnAmbyCode} \link{ReturnMrBayesAmbyCode}
#' @examples
#' ReturnMrBayesCode("A")

ReturnMrBayesCode <- function(base){
  if(base == "A")  return(1)
  if(base == "T")  return(2)
  if(base == "G")  return(3)
  if(base == "C")  return(4)
}

