phrynomics
==========

A collection of R code for dealing with SNP data and analyses. 

Check out the shiny app hosted on the UW Department of Statistics server
[shinyPhrynomics](https://rstudio.stat.washington.edu/shiny/phrynomics/)

##Installation

```r
devtools::install_github("bbanbury/phrynomics")
```


## Simple examples

Note: the library is still in active development and behaviour of the following
functions may well change in the future:

### Load a mock dataset


```r
library(phrynomics)
data(fakeData)
snpdata <- ReadSNP(fakeData)
```


![plot of chunk tree](http://i.imgur.com/YaxCqhe.png) 