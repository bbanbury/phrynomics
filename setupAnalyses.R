## Working file for setting up phylogenetic runs
##Barb Banbury  14 Nov 13
##See RAxMLanalysis for the tree analysis


################################################
########### Writing RAxML Data Files ###########
################################################

#Remove invariant sites and write data to file
source("~/phrynomics/trunk/phrynomicsFunctions.R")
setwd("~/Dropbox/phrynomics/Analyses/RAxMLruns")
snpFiles <- system("ls ~/Dropbox/phrynomics/c*", intern=T)
minSamples <- NULL
NuSNPs <- NULL
NuLoci <- NULL
for(i in sequence(length(snpFiles))) {
  filename <- strsplit(strsplit(snpFiles[i], "/")[[1]][6], ".", fixed=T)[[1]][1]
  snpList <- read.table(snpFiles[i], row.names=1, colClasses="character")
  snpList <- preProcess(snpList)
  loci <- dim(snpList)[2]
  initialLociLengths <- nchar(snpList[1,])
  splitdata <- SplitSNP(snpList)
  KeepVector <- apply(splitdata, 2, IsVariable)  #If True, then keep as variable
  breaks <- which(splitdata[1,] == " ")
  write.table(rbind(splitdata, KeepVector), file=paste("KeepVector", filename, ".txt", sep=""))

  newSNP <- cSNP(splitdata, KeepVector=KeepVector)
  newLoci <- length(which(newSNP[1,] != ""))  #number of new loci 

  SNPFile <- paste(dim(newSNP)[1], length(which(KeepVector == "TRUE")) - length(breaks))
  write(SNPFile, file=paste(filename, "noAmbigs.snps", sep=""))
  write.table(newSNP, file=paste(filename, "noAmbigs.snps", sep=""), append=TRUE, quote=FALSE, col.names=FALSE)  
  print(paste("wrote", i, "of", length(snpFiles)))
  minSamples <- c(minSamples, as.numeric(strsplit(filename, "\\D+", perl=T)[[1]][2]))
  NuSNPs <- c(NuSNPs, length(which(KeepVector == "TRUE")) - length(breaks))
  NuLoci <- c(NuLoci, newLoci)
}



#exploratory plots of data
plot(minSamples, log(NuSNPs), type="n")
text(minSamples, log(NuSNPs), labels=NuSNPs)
plot(minSamples, log(NuLoci), type="n")
text(minSamples, log(NuLoci), labels= NuLoci, col="red")



#Run RAxML on cluster
files <- system("ls c*", intern=T)
systemSeeds <- NULL
treeSeeds <- NULL
filenames <- NULL
models <- c("ASC_GTRGAMMA", "GTRGAMMA")
for(model in sequence(length(models))) {
  jobArgs <- NULL
  for(i in sequence(length(files))) {
    NumCores <- 16
    inputFile <- files[i]
    outgroup <- "GAWI1"
    systemSeed <- floor(runif(1, min=1, max=10^6))
    bootstrapReps <- 1000
    treeSeed <- floor(runif(1, min=1, max=10^6))
    outputName <- paste(gsub("\\D+$", "", files[i]), "_", models[model], sep="")
    systemCall <- paste("-T ", NumCores, " -s ", inputFile, " -f a -m ", models[model], " -o ", outgroup, " -x ", systemSeed, " -# ", bootstrapReps, " -p ", treeSeed, " -n ", outputName, sep="")
    jobArgs <- append(jobArgs, systemCall)

    #shopkeep
    systemSeeds <- c(systemSeeds, systemSeed) 
    treeSeeds <- c(treeSeeds, treeSeed) 
    filenames <- c(filenames, outputName) 
  }
  writeFileName <- paste(paste(strsplit(models[model], "")[[1]][1:3], collapse=""))
  save(systemSeeds, treeSeeds, filenames, file="seeds.rdata")
  write(jobArgs, file=paste(writeFileName, "_RAXML_JobArgs", sep=""))
}


#write shell script for each input file (run on imac)
setwd("~/Dropbox/UWstuff/phrynomics/Analyses/RAxMLruns")
files <- system("ls c*", intern=T)
models <- c("ASC_GTRCAT", "GTRCAT")
systemSeeds <- NULL
treeSeeds <- NULL
filenames <- NULL
setwd("~/Dropbox/UWstuff/phrynomics/Analyses/RAxMLruns2")
write(paste(""), file="raxmlSystemCalls.txt")
orderToGo <- c("c65p3", "c45p3", "c25p3", "c5p3", "c10p3", "c15p3", "c20p3", "c30p3", "c35p3", "c40p3", "c50p3", "c55p3", "c60p3", "c70p3")
for(model in sequence(length(models))){
  jobArg <- NULL
  for(i in sequence(length(snpFiles))){
    snpFile <- files[grep(orderToGo[i], files)]
    NumCores <- 8
    inputFile <- snpFile
    outgroup <- "GAWI1"
    systemSeed <- floor(runif(1, min=1, max=10^6))
    bootstrapReps <- 1000
    treeSeed <- floor(runif(1, min=1, max=10^6))
    outputName <- paste(gsub("\\D+$", "", snpFile), "_", models[model], sep="")
    filename <- NULL
    filename <- paste(models[model], "_", strsplit(strsplit(snpFile, "/")[[1]][9], "noAmbigs",  fixed=T)[[1]][1], sep="")
    jobArg <- paste("raxmlHPC-PTHREADS -T ", NumCores, " -s ../", inputFile, " -f a -m ", models[model], " -V -o ", outgroup, " -x ", systemSeed, " -# ", bootstrapReps, " -p ", treeSeed, " -n ", outputName, sep="")
    write(jobArg, file="raxmlSystemCalls.txt", append=TRUE)
    systemSeeds <- c(systemSeeds, systemSeed) 
    treeSeeds <- c(treeSeeds, treeSeed) 
    filenames <- c(filenames, outputName) 
  }
  save(systemSeeds, treeSeeds, filenames, file="seeds.rdata")
}

#to execute this new script, open in R
toRun <- scan(file="~/Dropbox/UWstuff/phrynomics/Analyses/RAxMLruns2/raxmlSystemCalls.txt", what="character", sep="\n")
for(i in sequence(length(toRun))){
  system(toRun[i])
}


################################################
###################### END #####################
################################################





################################################
########## Writing MrBayes Data Files ##########
################################################

setwd("~/Dropbox/UWstuff/phrynomics/Analyses/MrBayesRuns")
source("~/phrynomics/trunk/phrynomicsFunctions.R")
snpFiles <- system("ls ~/Dropbox/UWstuff/phrynomics/Analyses/RAxMLruns/c*", intern=T)
seedNumbers <- NULL
AllFilenames <- NULL
type <- c("variable", "all")
models <- c("ASC", "GTR")
for(model in sequence(length(models))){
  for(i in sequence(length(snpFiles))) {
    filename <- strsplit(strsplit(snpFiles[i], "/")[[1]][9], ".", fixed=T)[[1]][1]
    randomSeed <- floor(runif(1, min=1, max=100000))
    seedNumbers <- c(seedNumbers, randomSeed)
    snpList <- read.table(snpFiles[i], row.names=1, colClasses="character", skip=1)
    loci <- nchar(paste(snpList[1,], collapse=""))  #calculate number of loci BEFORE translating bases or the number of loci will be over estimated.  
    snpList2 <- TranslateBases(snpList, TRUE, "?")

    f <- paste(models[model], "_", filename, ".nex", sep="")  
    AllFilenames <- c(AllFilenames, f)
    write(paste("#NEXUS"), file=f)
    write(paste("[Written ", Sys.Date(), " via Rscript]", sep=""), file=f, append=TRUE)
    write(paste("BEGIN Data;"), file=f, append=TRUE)
    write(paste("   DIMENSIONS NTAX=", dim(snpList2)[1], " NCHAR=", loci, ";", sep=""), file=f, append=TRUE)
    write(paste("   FORMAT DATATYPE=Standard INTERLEAVE=no missing=?;", sep=""), file=f, append=TRUE)
    write(paste("Matrix"), file=f, append=TRUE)
    write.table(snpList2, file=f, append=TRUE, quote=FALSE, col.names=FALSE)  
    write(paste(""), file=f, append=TRUE)
    write(paste(";"), file=f, append=TRUE)
    write(paste("END;"), file=f, append=TRUE)
    write(paste(""), file=f, append=TRUE)
    write(paste("BEGIN mrbayes;"), file=f, append=TRUE)
    write(paste("outgroup GAWI1;"), file=f, append=TRUE)
    write(paste(""), file=f, append=TRUE)
    write(paste("set autoclose = yes;"), file=f, append=TRUE)
    write(paste("set seed = ", randomSeed, ";", sep=""), file=f, append=TRUE)
    write(paste("lset rates=gamma coding=", type[model], ";", sep=""), file=f, append=TRUE)
    write(paste("prset brlenspr=unconstrained:gammadir(1,1,1,1);"), file=f, append=TRUE)
    write(paste(""), file=f, append=TRUE)
    write(paste("mcmc nchains = 8 ngen = 20000000 printfreq = 200000 samplefreq = 10000;"), file=f, append=TRUE)
    write(paste("sumt burnin=500;", sep=""), file=f, append=TRUE)
    write(paste("sump burnin=500;", sep=""), file=f, append=TRUE)
    write(paste(""), file=f, append=TRUE)
    write(paste("END;"), file=f, append=TRUE)
    print(paste("wrote", i, "of", length(snpFiles)))
  }
  save(seedNumbers, AllFilenames, file="seeds.Rdata")
}


#write shell script for each input file
models <- c("ASC", "GTR")
snpFiles <- system("ls ~/Dropbox/phrynomics/Analyses/RAxMLruns/c*", intern=T)
for(model in sequence(length(models))){
  for(i in sequence(length(snpFiles))){
    filename <- NULL
    filename <- paste(models[model], "_", strsplit(strsplit(snpFiles[i], "/")[[1]][8], "noAmbigs",  fixed=T)[[1]][1], sep="")
    f <- paste("job", filename, ".sh", sep="")
    write(paste("#!/bin/sh"), file=f)
    write(paste(""), file=f, append=TRUE)
    write(paste('#PBS -N "', filename, '"', sep=""), file=f, append=TRUE)
    write(paste("#PBS -l nodes=1:ppn=12,feature=12core,mem=22gb,walltime=500:00:00"), file=f, append=TRUE)
    write(paste("#PBS -o /gscratch/leache/Barb/jobLogs/"), file=f, append=TRUE)
    write(paste("#PBS -d /gscratch/leache/Barb/asc_MrBayesjobs/", models[model], "/", strsplit(strsplit(snpFiles[i], "/")[[1]][8], "noAmbigs",  fixed=T)[[1]][1], sep=""), file=f, append=TRUE)
    write(paste("module load gcc_4.1.2-ompi_1.6.2"), file=f, append=TRUE)
    write(paste("mpirun -np 12 /gscratch/leache/Barb/mrbayes_3.2.2/src/mb ", filename, "noAmbigs.nex ", "> log", filename, ".txt", sep=""), file=f, append=TRUE)
  }
}

################################################
###################### END #####################
################################################





################################################
########### Writing SNAPP Data Files ###########
################################################

#file <- ("~/Dropbox/UWstuff/phrynomics/Analyses/RAxMLruns/c65p3noAmbigs.snps")
# file <- ("~/Dropbox/UWstuff/phrynomics/Analyses/RAxMLruns/c5p3noAmbigs.snps")  #old file that makes 36 final SNPs
file <- ("~/Dropbox/UWstuff/phrynomics/Analyses/SNAPP/c20p25.unlinked_snps")

source("~/phrynomics/trunk/phrynomicsFunctions.R")

Data1 <- read.table(file, row.names=1, colClasses="character", skip=1)
#Data2 <- removeOutgroups(Data1) #all removed for new unlinked file
Data2 <- Data1
#Data3 <- RemoveMissingSpeciesLoci(SplitSNP(Data2))
Data3 <- Data2
Data4 <- RemoveNonBinary(Data3)
#Data5 <- TakeSingleSNPfromEachLocus(Data4)[[1]]
Data5 <- Data4
Data6 <- TranslateBases(Data5, translateMissing=FALSE, ordered=TRUE)

setwd("~/Dropbox/UWstuff/phrynomics/Analyses/SNAPP")
SNAPPfile <- paste(dim(Data6)[1], nchar(Data6[1,]))
# write(SNAPPfile, file="allPhrynosomaSNAPP.noStep3.snps")
# write.table(Data6, file="allPhrynosomaSNAPP.noStep3.snps", append=TRUE, quote=FALSE, col.names=FALSE)  
write(SNAPPfile, file="translated.c20p25.unlinked_snps")
write.table(Data6, file="translated.c20p25.unlinked_snps", append=TRUE, quote=FALSE, col.names=FALSE)  

#remove individuals with missing data (not just species, ut now individuals)
indToRemove <- NULL
for(i in sequence(dim(Data6)[1])){
  if(Data6[i,])
}




################################################
###################### END #####################
################################################




################################################
######  Comparison using invariant Sites  ######
################################################

#checking allSites using sims
snpList <- read.table("sim.phy", row.names=1, colClasses="character", skip=1)
  loci <- dim(snpList)[2]
  initialLociLengths <- nchar(snpList[1,])
  splitdata <- SplitSNP(snpList)
  KeepVector <- apply(splitdata, 2, IsVariable)  #If True, then keep as variable
  breaks <- which(splitdata[1,] == " ")
  newSNP <- cSNP(splitdata, KeepVector=KeepVector)
  newLoci <- length(which(newSNP[1,] != ""))  #number of new loci 
  SNPFile <- paste(dim(newSNP)[1], length(which(KeepVector == "TRUE")) - length(breaks))
  write(SNPFile, file="red.sim.phy")
  write.table(newSNP, file="red.sim.phy", append=TRUE, quote=FALSE, col.names=FALSE)  

#checking allSites using c55 dataset and made up invariant sites
snpList <- read.table("~/Dropbox/UWstuff/phrynomics/Analyses/RAxMLruns2/c55p3noAmbigs.snps", row.names=1, colClasses="character", skip=1)
extra1 <- as.character(paste(rep("A", 500), collapse=""))
extra2 <- as.character(paste(rep("G", 500), collapse=""))
extra3 <- as.character(paste(rep("T", 500), collapse=""))
extra4 <- as.character(paste(rep("C", 500), collapse=""))
snpList <- cbind(snpList, extra1, extra2, extra3, extra4, stringsAsFactors=FALSE)
loci <- dim(snpList)[2]
initialLociLengths <- nchar(snpList[1,])
SNPFile <- paste(dim(snpList)[1], sum(initialLociLengths))
write(SNPFile, file="c55.extraData2")
write.table(snpList, file="c55.extraData2", append=TRUE, quote=FALSE, col.names=FALSE)  


#after simulating data along a random 73 taxon tree
#get ratio to the same variant/total sites as c45 dataset
c45data <- read.table("~/Dropbox/UWstuff/phrynomics/c45p3.phy", colClasses="character", row.names=1, skip=1)
c45pruned <- read.table("~/Dropbox/UWstuff/phrynomics/Analyses/RAxMLruns2/c45p3noAmbigs.snps", colClasses="character", row.names=1, skip=1)
ratio <- sum(nchar(c45pruned[1,])) / nchar(c45data[1,])   #16.6%

library(phangorn)
tree <- rtree(73)
tree <- root(tree, "t10")
tree$edge.length <- tree$edge.length/400
data <- simSeq(tree, l=10000, type="DNA", bf=c(.25,.25,.25,.25), Q=rep(1,6))
data <- toupper(as.character(data))
initialLociLengths <- length(data[1,])
data2 <- cSNP(data)
SNPFile <- paste(dim(data2)[1], initialLociLengths)
write(SNPFile, file="c45.sim")
write.table(data2, file="c45.sim", append=TRUE, quote=FALSE, col.names=FALSE)  

KeepVector <- apply(data, 2, IsVariable)  #If True, then keep as variable
newSNP <- cSNP(data, KeepVector=KeepVector)
newLociLength <- nchar(newSNP[1,])
newLociLength/initialLociLengths  #16.3%
SNPFile <- paste(dim(newSNP)[1], newLociLength)
write(SNPFile, file="c45.red.sim")
write.table(newSNP, file="c45.red.sim", append=TRUE, quote=FALSE, col.names=FALSE)  

## setup three RAxML runs on two sim files (full, nonASC, ASC)
setwd("~/Dropbox/UWstuff/phrynomics/Analyses/RAxML.allSites")
fullRun <- paste("raxmlHPC-PTHREADS -T 8 -s c45.sim -f a -m GTRCAT -o t10 -x 320234 -# 4 -p 30947 -n c45.full")
igfullRun <- paste("raxmlHPC-PTHREADS -T 8 -s c45.sim -f a -m GTRGAMMAI -o t10 -x 2439085 -# 4 -p 981437 -n c45.igfull")
gfullRun <- paste("raxmlHPC-PTHREADS -T 8 -s c45.sim -f a -m GTRGAMMA -o t10 -x 249085 -# 4 -p 98437 -n c45.gfull")
nonASCrun <- paste("raxmlHPC-PTHREADS -T 8 -s c45.red.sim -f a -m GTRCAT -o t10 -x 64831 -# 4 -p 186324 -n c45.nonASC")
ASCrun <- paste("raxmlHPC-PTHREADS -T 8 -s c45.red.sim -f a -m ASC_GTRCAT -V -o t10 -x 192873 -# 4 -p 75984 -n c45.ASC")

allRuns <- c(fullRun, igfullRun, gfullRun, nonASCrun, ASCrun)
for(i in sequence(length(allRuns))){
  system(allRuns[i], intern=FALSE)
}
























