##  --------------------------------  ##
##  ---      shinyPhrynomics     ---  ##
##  ---    written by B.Banbury  ---  ##
##  ---     v1.3  8August2014    ---  ##
##  --------------------------------  ##


##  ---         Load Code        ---  ##
library(shiny)
library(ape)
library(phangorn)
library(devtools)
#devtools::install_github("bbanbury/phrynomics")
#library(phrynomics)
load_all("~/phrynomics/branches/OOP")
library(shinyIncubator)  #maybe use for the progress bar

vers <- "v1.3"


##  ---     Server Functions     ---  ##


convertNexDataToPhyData <- function(nexData) {
  phyData <- data.frame(matrix(lapply(lapply(nexData, toupper), paste, collapse="")))
  rownames(phyData) <- names(nexData)
  return(phyData)
}



##  ---   Server Communication   ---  ##
shinyServer(function(input, output) {
  
  # For files that are over 5MB, we can set the acceptance to 15. 
  options(shiny.maxRequestSize=15*1024^2)

  output$phrynoversion <- renderText({
    paste("Made with Phrynoversion", vers)
  })


  ##  ---   Initial Data   ---  ##

  initializeTable <- reactive({
    if(is.null(input$OrigData))
      return(NULL)
    inFile <- input$OrigData
    fileSize <- inFile$size
    inputFileType <- FileFormat(inFile$datapath)
    if(is.null(inputFileType))
      return(NULL)
    if(inputFileType == "phy")
      initializeTable <- ReadSNP(inFile$datapath)
    if(inputFileType == "nex") 
      initializeTable <- convertNexDataToPhyData(read.nexus.data(inFile$datapath))  #change this over to ReadSNP too
    return(initializeTable)
  })
  
  ##  ---   Modified Data   ---  ##

  contents <- reactive({
    textOutput <- NULL
    results <- initializeTable()$data
    if (is.null(results))
      return(NULL)
    #if(input$minInd > 1){
#      results <- ReduceMinInd(results, "loci", threshold=input$minInd)
    #}
    if(input$rmInvSites)
      results <- RemoveInvariantSites(results)
    if(input$rmNonBin)
      results <- RemoveNonBinary(results)
    if(input$takeRandom)
      results <- TakeSingleSNPfromEachLocus(results)[[1]]
    if(input$transSNAPP){
      if(all(apply(SplitSNP(results), 2, IsBinary), apply(SplitSNP(results), 2, IsVariable))){
        textOutput <- "Calculating..."
        results <- TranslateBases(results, translateMissingChar=input$transformChar, ordered=TRUE)
        textOutput <- "Transformed to SNAPP"
      }
      else
        textOutput <- c(textOutput, "Some sites are not binary or variable, can not translate yet!")        
    }
    if(input$transMrBayes) {
      results <- TranslateBases(results, translateMissingChar=input$transformChar, ordered=FALSE)
      textOutput <- c(textOutput, "Preview is not supported for MrBayes transformations")
      textOutput <- c(textOutput, "Download results should still work! Check it!")
    }
  return(list(ReadSNP(results), textOutput))    
  })

  output$slider <- renderUI({
    maxTax <- contents()[[1]]$ntax
    if(is.null(maxTax))
      maxTax <- 10
    sliderInput("minInd", "Minimum number of individuals with sequence data required for a locus to be included in the data set.", min=1, max=maxTax, value=(0.25*maxTax), round=TRUE, ticks=FALSE)
  })

  ##  ---    Summary Tab    ---  ##

  output$summaryData <- renderUI({
    if(is.null(contents()[[1]]))
      return(NULL)
    numLoci <- contents()[[1]]$nloci
    space <- ""
    a <- paste("Number of taxa:", contents()[[1]]$ntax)
    b <- paste("Number of loci:", numLoci)
    if(numLoci == 1)
      b <- paste("One concatenated dataset")
    c <- paste("Number of sites:", sum(contents()[[1]]$nsites))
    d <- paste("Average number of sites per locus:", round(mean(contents()[[1]]$nsites), digits=3))
    e <- paste("Minimum number of sites per locus:", min(contents()[[1]]$nsites))
    f <- paste("Maximum number of sites per locus:", max(contents()[[1]]$nsites))
    HTML(paste(a, b, c, space, d, e, f, sep="<br/>"))
  })


  ##  ---   Main Data Tab   ---  ##

  output$communicationWindow <- renderText({
    if(is.null(contents()[[2]]))
      return(NULL)
    contents()[[2]]
  })

  output$OrigDataStats <- renderText({
    forStats <- initializeTable()
    if(is.null(forStats))
      return(NULL)
    return(paste("Original Data contains", forStats$ntax, "taxa,", forStats$nloci, "loci, and", sum(forStats$nsites), "sites"))
  })

  output$resultsDataStats <- renderText({
    resultsStats <- contents()[[1]]
    if(is.null(resultsStats))
      return(NULL)
    return(paste("Transformed Data contains", resultsStats$ntax, "taxa,", resultsStats$nloci, "loci, and", sum(resultsStats$nsites), "sites"))
  })


 output$inputPreview <- renderTable({  
    prev <- head(initializeTable()$data, n=input$obs)
    newDim <- min(sum(nchar(prev[1,])), input$snpobs) 
    if(!is.null(prev))
      return(cSNP(SplitSNP(prev)[,1:newDim]))
  }, include.colnames=FALSE, size=1)
  
  output$resultsPreview <- renderTable({   
    outPrev <- head(contents()[[1]]$data, n=input$obs)
    newDim2 <- min(sum(nchar(outPrev[1,])), input$snpobs)      
    if(input$transMrBayes)
      return(NULL)
    if(!is.null(outPrev))
      outPrev <- cSNP(SplitSNP(outPrev)[,1:newDim2])
  }, include.colnames=FALSE, size=1)
  
  output$downloadPhy <- downloadHandler(
    filename = function() {paste(input$datasetName, ".txt", sep="")},
    content = function(file) {WriteSNP(contents()[[1]], file, format="phylip")})

  output$downloadNex <- downloadHandler(
    filename = function() { paste(input$datasetName, ".nex", sep="") },
    content = function(file) {WriteSNP(contents()[[1]], file, format="nexus", missing=input$transformChar)})


  ##  ---   Plots Data Tab   ---  ##

  output$summaryPlot1 <- renderPlot({
    if(is.null(contents()[[1]]$data))
      return(NULL)
    x <- contents()[[1]]
    plot(x$nsites, type="n", ylab="Number Sites", xlab="Locus", xaxt="n")
    text(1:length(x$nsites), x$nsites, labels=1:length(x$nsites), col="blue")
    segments(0, mean(x$nsites), length(x$nsites)+1, mean(x$nsites), col="red")
    text(0.1+(0.1*length(x$nsites)), mean(x$nsites)+(.02*max(x$nsites)), "mean", col="red")
    title(main="Number of Sites Per Locus")
  })

  output$summaryPlot2 <- renderPlot({
    if(is.null(contents()[[1]]$data))
      return(NULL)
    x <- contents()[[1]]
    plot(density(x$nsites), xlab="Number Sites", ylab="Frequency", main="", col="blue")
    segments(mean(x$nsites), -0.1, mean(x$nsites), max(1, max(density(x$nsites)$y)), col="red")
    title(main="Frequency of Sites per Locus")
  })




})

