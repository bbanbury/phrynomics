summary.snp <- function (x, ...){
  cat("\nSNP dataset:", deparse(substitute(x)), "\n")
  cat(paste("\nSNP dataset contains",  x$ntax, "taxa,", x$nloci, "loci, and", sum(x$nsnps), "SNPs\n\n"))
  cat("\tNumber of taxa:", x$ntax, "\n")
  cat("\tNumber of loci:", x$nloci, "\n")
  cat("\tNumber of snps:", sum(x$nsnps), "\n")
  cat("\n\tAverage number of SNPs per locus:", mean(x$nsnps), "\n")
  cat("\tMinimum number of SNPs per locus:", min(x$nsnps), "\n")
  cat("\tMaximum number of SNPs per locus:", max(x$nsnps), "\n")
}