print.snp <- function (x, ...){
  cat(paste("\nSNP dataset contains",  x$ntax, "taxa,", x$nloci, "loci, and", sum(x$nsnps), "SNPs\n\n"))
  cat("Taxa include:\n")
  cat(paste("\t", paste(rownames(x$data)[1:min(6, ntax)], collapse = ", "), ", ...\n", sep = ""))
}
