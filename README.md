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


##Examples

Phrynomics can do a lot of dataset manipulations and simple calculations. 

Load a mock dataset. Phrynomics accepts either phylip or nexus sequential (ie, not interleaved) formatted files.  

```r
library(phrynomics)
data(fakeData)
snpdata <- ReadSNP(fakeData)
```

Phrynomics loads datasets into the class "snp". Typing the object into the console will only return a summary of the data, however you can see the data itself using the $ operator. 

There are lots of different measurements (calculate missing data, allele frequencies, number of sites per locus, etc) and also lots of ways to manipulate the data (translating bases, taking a single random SNP from each locus, removing invariant or nonbinary sites, remove individuals, add a flag to species names, etc) and ways to visualize the data in tables or plots (plot, making a heatmap of missing data, making a plot of missing data per minimum individuals, minimum individuals tables, etc.). 

### Some simple measures of your dataset:
```r
CalculateMissingData(snpdata, "loci")  
DataOverlap(snpdata) 
GetNumberOfSitesForLocus(snpdata) 
GetSpecies(rownames(snpdata$data)) 
MakeMinIndTable(snpdata, showEvery=1)
```

There are other simple measures that you can do per site (or use apply to do a whole dataset).  
```r
splits <- SplitSNP(snpdata$data)
apply(splits, 2, GetBaseFrequencies)
apply(splits, 2, ReturnUniqueBases)
apply(splits, 2, IsBinary)
apply(splits, 2, IsMissing)
apply(splits, 2, IsVariable)
```

### Examples of data transformation
```r
# Remove loci/sites with too much missing data
ReduceMinInd(snpdata, threshold=0.9)

# Remove Invariant/nonbinary sites or loci/sites with too much missing data
RemoveInvariantSites(snpdata)
RemoveNonBinary(snpdata)

# Add a flag to taxa to unify them as a population
AddAFlag(snpdata, flagToAdd="north", taxa=c("in1", "in1"))

# Convert Missing data
ConvertMissingData(snpdata)

#Take a single random site from each locus
TakeSingleSNPfromEachLocus(snpdata)

#Translate bases from alpha to numeric (used for MrBayes Mkv and SNAPP)
TranslateBases(snpdata)
```

###Plotting functions
```r
plot(snpdata)
plotHeatmap(snpdata)
plotMissing(snpdata)
```

###Exporting Data

Data can be exported in various formats (phylip or nexus) and for various phylogenetic programs.  

```r
# Write nexus or phylip formatted files
WriteSNP(snpdata, file="mydata.phy")
WriteSNP(snpdata, file="mydata.nex", format="nexus")
```

###Example Workflows 

Here are a few very minimal examples of what you can do with phrynomics from start to finish.  

```r
# Translate to MrBayes MkV formatting
snps <- ReadSNP("MyFileOfSNPs.phy")
snps <- RemoveInvariantSites(snps)
snps <- TranslateBases(snps, translateMissing=FALSE, ordered=FALSE)
WriteSNP(snps, file="snps.nex", format="nexus")

# Prepare data for RAxML
snps <- ReadSNP("MyFileOfSNPs.phy")
snps <- RemoveInvariantSites(snps)
snps <- RemoveNonBinary(snps)
WriteSNP(snps, file="snps.phy", format="phylip")

# Prepare data for SNPAPP
snps <- ReadSNP("MyFileOfSNPs.phy")
snps <- RemoveNonBinary(snps)
snps <- TakeSingleSNPfromEachLocus(snps)$snpdata
snps <- TranslateBases(snps, translateMissingChar="?", ordered=TRUE)
WriteSNP(snps, file="snps.nex", format="nexus", missing="?")

# Prepare data for TreeMix
snps <- ReadSNP("MyFileOfSNPs.phy")
snps <- RemoveNonBinary(RemoveInvariantSites(snps))
flagged.snps <- AddAFlag(snps, flagToAdd="WEST", taxa=c("in1", "in2", "in3"))
flagged.snps <- AddAFlag(flagged.snps, flagToAdd="EAST", taxa=c("in4", "in5", "in6"))
ExportPops(flagged.snps, subsets=c("WEST", "EAST"))
```

