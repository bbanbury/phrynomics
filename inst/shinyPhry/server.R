##  --------------------------------  ##
##  ---      shinyPhrynomics     ---  ##
##  ---    written by B.Banbury  ---  ##
##  ---     v1.2  10June2014    ---  ##
##  --------------------------------  ##


##  ---         Load Code        ---  ##
library(shiny)
library(ape)
library(phangorn)
library(devtools)
devtools::install_github("bbanbury/phrynomics")
library(phrynomics)

vers <- "v1.2"


##  ---     Server Functions     ---  ##
fileFormat <- function(file){
  # returns "phy" or "nex" depending on the file type
  format <- NULL
  if(any(grep("Error", try(read.table(file, skip=GetLinesToSkip(file)+1), silent=TRUE), ignore.case=TRUE))) {
    if(!any(grep("Error", try(read.nexus.data(file)), ignore.case=TRUE)))
      format <- "nex"
  } else{
    format <- "phy"
  }
  return(format)
}

convertNexDataToPhyData <- function(nexData) {
  phyData <- data.frame(matrix(lapply(lapply(nexData, toupper), paste, collapse="")))
  rownames(phyData) <- names(nexData)
  return(phyData)
}

convertPhyDataToNexData <- function(phyData) {
  phyData <- SplitSNP(phyData)
  if(any(phyData[1,] == " "))
    phyData <- SplitSNP(phyData)[,-which(SplitSNP(phyData)[1,] == " ")]
  nexData <- vector("list", dim(phyData)[1])
  for(i in 1:length(nexData)){
    nexData[[i]] <- phyData[i,] #nexus parser doesn't like spaces between loci
    names(nexData)[i] <- rownames(phyData)[i]
  }
  return(nexData)
}

GetNumberSNPs <- function(taxon){
  return(nchar(paste(taxon, collapse="")))
}

WriteNexus <- function(phyData, file, missing) {
  if(class(phyData) == "data.frame" || class(phyData) == "matrix")
    if(dim(phyData)[2] > 1)
      phyData <- apply(phyData, 1, paste, collapse="")
  nchars <- min(sapply(phyData, GetNumberSNPs))
  write(paste("#NEXUS"), file)
  write(paste("[Written ", Sys.Date(), " via shinyPhrynomics]", sep=""), file, append=TRUE)
  write(paste("BEGIN Data;"), file, append=TRUE)
  write(paste("   DIMENSIONS NTAX=", length(phyData), " NCHAR=", nchars, ";", sep=""), file, append=TRUE)
  write(paste("   FORMAT DATATYPE=Standard INTERLEAVE=no missing=", missing, ";", sep=""), file, append=TRUE)
  write(paste("Matrix"), file, append=TRUE)
  write.table(phyData, file, append=TRUE, quote=FALSE, col.names=FALSE)  
  write(paste(""), file, append=TRUE)
  write(paste(";"), file, append=TRUE)
  write(paste("END;"), file, append=TRUE)
  write(paste(""), file, append=TRUE)
}

sum_along <- function(x){
  sums <- rep(NA, length(x))
  for(i in sequence(length(sums))){
    sums[i] <- sum(x[1:i])
  }
return(sums)
}

##  ---   Server Communication   ---  ##
shinyServer(function(input, output) {
  output$phrynoversion <- renderText({
    paste("Made with Phrynoversion", vers)
  })

  initializeTable <- reactive({
    inFile <- input$SNPdataset
    inputFileType <- fileFormat(inFile$datapath)
    if (is.null(inputFileType))
      return(NULL)
    if(inputFileType == "phy")
      initializeTable <- read.table(inFile$datapath, row.names=1, skip=1, stringsAsFactors=FALSE)
    if(inputFileType == "nex") 
      initializeTable <- convertNexDataToPhyData(read.nexus.data(inFile$datapath))
    return(initializeTable)
  })
  
  contents <- reactive({
    textOutput <- NULL
    results <- initializeTable()
    if (is.null(results))
      return(NULL)
    if(input$rmInvSites)
      results <- RemoveInvariantSites(results)
    if(input$rmNonBin)
      results <- RemoveNonBinary(results)
    if(input$takeRandom)
      results <- TakeSingleSNPfromEachLocus(results)[[1]]
    if(input$transSNAPP){
      results <- TranslateBases(results, translateMissingChar=input$transformChar, ordered=TRUE)
      textOutput <- c(textOutput, "Transformed to SNAPP")
    }
    if(input$transMrBayes) {
      results <- TranslateBases(results, translateMissingChar=input$transformChar, ordered=FALSE)
      textOutput <- c(textOutput, "Preview is not supported for MrBayes transformations")
      textOutput <- c(textOutput, "Download results should still work! Check it!")
    }
  return(list(results, textOutput))    
  })

  output$OrigDataStats <- renderText ({
    forStats <- initializeTable()
    if(is.null(forStats))
      return(NULL)
    return(paste("Original Data contains", dim(forStats)[1], "taxa,", dim(forStats)[2], "sites, and", GetNumberSNPs(forStats[1,]), "SNPs"))
  })

  output$resultsDataStats <- renderText ({
    resultsStats <- contents()[[1]]
    if(is.null(resultsStats))
      return(NULL)
    return(paste("Transformed Data contains", dim(resultsStats)[1], "taxa,", dim(resultsStats)[2], "sites, and", GetNumberSNPs(resultsStats[1,]), "SNPs"))
  })

  output$inputPreview <- renderTable({  
    prev <- head(initializeTable(), n=input$obs)
    newDim <- min(sum(nchar(prev[1,])), input$snpobs) 
    if(!is.null(prev))
      return(cSNP(SplitSNP(prev)[,1:newDim]))
  }, include.colnames=FALSE, size=1)
  
  output$resultsPreview <- renderTable({   
    outPrev <- head(contents()[[1]], n=input$obs)
    newDim2 <- min(sum(nchar(outPrev[1,])), input$snpobs)      
    if(input$transMrBayes)
      return(NULL)
    if(!is.null(outPrev))
      outPrev <- cSNP(SplitSNP(outPrev)[,1:newDim2])
  }, include.colnames=FALSE, size=1)
  
  output$downloadPhy <- downloadHandler(
    filename = function() { paste(input$datasetName, ".txt", sep="") },
    content = function(file) {
      fileContent <- paste(dim(contents()[[1]])[1], sum(nchar(contents()[[1]][1,])))
      write(fileContent, file)
      write.table(contents()[[1]], file, quote=FALSE, append=TRUE, col.names=FALSE)
    })

  output$downloadNex <- downloadHandler(
    filename = function() { paste(input$datasetName, ".nex", sep="") },
    content = function(file) {
      fileContent <- WriteNexus(contents()[[1]], file, missing=input$transformChar)
  })

output$communicationWindow <- renderText ({
  if(is.null(contents()[[2]]))
    return(NULL)
  contents()[[2]]
  })

})

