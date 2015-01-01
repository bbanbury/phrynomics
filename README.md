phrynomics
==========
![](http://barbbanbury.info/barbbanbury/Research_Projects_files/phrynoHead.jpg) 

A collection of R code for dealing with SNP data and analyses. 

Check out the shiny app hosted on the UW Department of Statistics server
[shinyPhrynomics](https://rstudio.stat.washington.edu/shiny/phrynomics/)

##Installation

```r
devtools::install_github("bbanbury/phrynomics")
```


## Phrynomics can do a lot of dataset manipulations and simple calculations. 

### Load a mock dataset. Phrynomics accepts either phylip or nexus sequential (ie, not interleaved) formatted files.  

```r
library(phrynomics)
data(fakeData)
snpdata <- ReadSNP(fakeData)
```

### Phrynomics loads datasets into the class "snp". There are lots of different measurements (calculate missing data, allele frequencies, number of sites per locus, etc) and also lots of ways to manipulate the data (translating bases, taking a single random SNP from each locus, removing invariant or nonbinary sites, remove whole species (multiple individuals), add a flag to species names, etc) and ways to visualize the data in tables or plots (plot, making a heatmap of missing data, making a plot of missing data per minimum individuals, minimum individuals tables, etc.). 

### Examples of simple measures of your dataset:
```r
CalculateMissingData(snpdata, "loci")
DataOverlap(snpdata)
```

