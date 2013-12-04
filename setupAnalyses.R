## Working file for setting up phylogenetic runs
##Barb Banbury  14 Nov 13
##See RAxMLanalysis for the tree analysis


################################################
########### Writing RAxML Data Files ###########
################################################

#Remove invariant sites and write data to file
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

################################################
###################### END #####################
################################################





################################################
########## Writing MrBayes Data Files ##########
################################################

setwd("~/Dropbox/UWstuff/phrynomics/Analyses/MrBayesRuns")
source('~/Dropbox/UWstuff/phrynomics/Analyses/Rcode/PhrynoNucCodes.R', chdir = TRUE)
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



























