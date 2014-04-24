##  --------------------------------  ##
##  ---      shinyPhrynomics     ---  ##
##  ---    written by B.Banbury  ---  ##
##  ---     v1.0  22April2014    ---  ##
##  --------------------------------  ##


##  ---         Load Code        ---  ##
library(shiny)
vers <- "v1.0"
loadSourceFiles <- function(path){
  filesToLoad <- list.files(path)
  sapply(paste(path, filesToLoad, sep=""), source)
}
loadSourceFiles("../../R/")

##  ---     Server Functions     ---  ##
fileFormat <- function(file){
  # returns "phy" or "nex" depending on the file type
  format <- NULL
  if(any(grep("Error", try(read.table(file, skip=1), silent=TRUE), ignore.case=TRUE))) {
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
  nchars <- min(apply(phyData, 1, GetNumberSNPs))
  write(paste("#NEXUS"), file)
  write(paste("[Written ", Sys.Date(), " via shinyPhrynomics]", sep=""), file, append=TRUE)
  write(paste("BEGIN Data;"), file, append=TRUE)
  write(paste("   DIMENSIONS NTAX=", dim(phyData)[1], " NCHAR=", nchars, ";", sep=""), file, append=TRUE)
  write(paste("   FORMAT DATATYPE=Standard INTERLEAVE=no missing=", missing, ";", sep=""), file, append=TRUE)
  write(paste("Matrix"), file, append=TRUE)
  write.table(phyData, file, append=TRUE, quote=FALSE, col.names=FALSE)  
  write(paste(""), file, append=TRUE)
  write(paste(";"), file, append=TRUE)
  write(paste("END;"), file, append=TRUE)
  write(paste(""), file, append=TRUE)
}


##  ---   Server Communication   ---  ##
shinyServer(function(input, output) {
  output$phrynoversion <- renderText({
    paste("Made with Phrynoversion", vers)
  })

  initializeTable <- reactive({
    inFile <- input$SNPdataset
    inputFileType <- fileFormat(inFile$datapath)
    if (is.null(inFile))
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
      
  output$inputPreview <- renderTable({  
    head(initializeTable(), n=input$obs)
  }, include.colnames=FALSE, size=1)
  
  output$resultsPreview <- renderTable({   
    outPrev <- head(contents()[[1]], n=input$obs)
#    maxnchars <- max(apply(outPrev, 1, GetNumberSNPs))
#    spacesToAdd <- maxnchars - apply(outPrev, 1, GetNumberSNPs)
#    for(i in 1:dim(outPrev)[1]){
#      spacesToAdd[i] <- paste(rep("0", spacesToAdd[i]), collapse="")
#      #outPrev[i] <- paste(outPrev[i,], spacesToAdd[i], sep="")
#    }
#    save(outPrev, file="outPrev.rdata")
    if(input$transMrBayes)
      return(NULL)
    if(!is.null(outPrev))
      outPrev <- cSNP(SplitSNP(outPrev))
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

