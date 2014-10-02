##  --------------------------------  ##
##  ---      shinyPhrynomics     ---  ##
##  ---    written by B.Banbury  ---  ##
##  ---     v1.3  15August2014    ---  ##
##  --------------------------------  ##


##  ---         Load Code        ---  ##
library(shiny)
library(ape)
library(phangorn)
library(devtools)
#devtools::install_github("bbanbury/phrynomics")
library(phrynomics)
#load_all("~/phrynomics/branches/OOP")

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
    if(initializeTable$nloci == "1")
      initializeTable <- ReadSNP(SplitSNP(l$data))
    return(initializeTable)
  })
  
  ##  ---   Modified Data   ---  ##

  contents <- reactive({
    textOutput <- NULL
    results <- initializeTable()$data
    if (is.null(results))
      return(NULL)
    if(input$minInds != ""){
      results <- ReduceMinInd(results, "loci", threshold=input$minInds)
    }
    results <- RemoveInvariantSites(results)
  return(ReadSNP(results))
  })

  myplots <- reactive({
    if(is.null(contents()$data))
      return(NULL)
    x <- contents()
    if(input$plotIn == ""){
      return(NULL)
    }
    if(input$plotIn == "Per Locus"){
      plot(x$nsites, type="n", ylab="Number Sites", xlab="Locus", xaxt="n")
      text(1:length(x$nsites), x$nsites, labels=1:length(x$nsites), col="blue")
      segments(0, mean(x$nsites), length(x$nsites)+1, mean(x$nsites), col="red")
      text(0.1+(0.1*length(x$nsites)), mean(x$nsites)+(.02*max(x$nsites)), "mean", col="red")
      title(main="Number of Sites Per Locus")
    }
    if(input$plotIn == "Frequency"){
      plot(density(x$nsites), xlab="Number Sites", ylab="Frequency", main="", col="blue")
      segments(mean(x$nsites), -0.1, mean(x$nsites), max(1, max(density(x$nsites)$y)), col="red")
      title(main="Frequency of Sites per Locus")
    }
    if(input$plotIn == "Missing Data"){
      plotMissing(x)
      #segments(mean(x$nsites), -0.1, mean(x$nsites), max(1, max(density(x$nsites)$y)), col="red")
      #title(main="Frequency of Sites per Locus")
    }
    if(input$plotIn == "Heatmap"){
      plotHeatmap(x)
      #segments(mean(x$nsites), -0.1, mean(x$nsites), max(1, max(density(x$nsites)$y)), col="red")
      #title(main="Frequency of Sites per Locus")
    }
  })

  myplots2 <- function(){
    if(is.null(contents()$data))
      return(NULL)
    x <- contents()
    if(input$plotIn == ""){
      return(NULL)
    }
    if(input$plotIn == "Per Locus"){
      plot(x$nsites, type="n", ylab="Number Sites", xlab="Locus", xaxt="n")
      text(1:length(x$nsites), x$nsites, labels=1:length(x$nsites), col="blue")
      segments(0, mean(x$nsites), length(x$nsites)+1, mean(x$nsites), col="red")
      text(0.1+(0.1*length(x$nsites)), mean(x$nsites)+(.02*max(x$nsites)), "mean", col="red")
      title(main="Number of Sites Per Locus")
    }
    if(input$plotIn == "Frequency"){
      plot(density(x$nsites), xlab="Number Sites", ylab="Frequency", main="", col="blue")
      segments(mean(x$nsites), -0.1, mean(x$nsites), max(1, max(density(x$nsites)$y)), col="red")
      title(main="Frequency of Sites per Locus")
    }
    if(input$plotIn == "Missing Data"){
      plotMissing(x)
      #segments(mean(x$nsites), -0.1, mean(x$nsites), max(1, max(density(x$nsites)$y)), col="red")
      #title(main="Frequency of Sites per Locus")
    }
    if(input$plotIn == "Heatmap"){
      plotHeatmap(x)
      #segments(mean(x$nsites), -0.1, mean(x$nsites), max(1, max(density(x$nsites)$y)), col="red")
      #title(main="Frequency of Sites per Locus")
    }
  }




  ##  ---   Rendering Output   ---  ##


  
  output$resultsPreview <- renderTable({   
    outPrev <- head(contents()$data, n=input$obs)
    newDim2 <- min(sum(nchar(outPrev[1,])), input$snpobs)      
    if(input$transMrBayes)
      return(NULL)
    if(!is.null(outPrev))
      outPrev <- cSNP(SplitSNP(outPrev)[,1:newDim2])
  }, include.colnames=FALSE, size=1)
  

  output$dRAxML <- downloadHandler(
    filename = function() {paste(input$datasetName, ".RAxML.phy", sep="")},
    content = function(file) {
      results <- RemoveNonBinary(contents())
      WriteSNP(results, file, format="phylip")})

  output$dMrBayes <- downloadHandler(
    filename = function() { paste(input$datasetName, ".MrBayes.nex", sep="") },
    content = function(file) {
      results <- TranslateBases(contents(), translateMissingChar="?", ordered=FALSE)     
      WriteSNP(results, file, format="nexus", missing="?")
    })

  output$dsnapp <- downloadHandler(
    filename = function() { paste(input$datasetName, ".SNAPP.nex", sep="") },
    content = function(file) {
      results <- RemoveNonBinary(contents())
      results <- TakeSingleSNPfromEachLocus(results)[[1]]
      results <- TranslateBases(results, translateMissingChar="?", ordered=TRUE)
      WriteSNP(results, file, format="nexus", missing="?")
    })


  ##  ---   Plots Data Tab   ---  ##

  output$summaryPlot1 <- renderPlot({
    myplots()
  })

  output$downloadSummPlot <- downloadHandler(
    filename = function() {paste0(input$datasetName, "summary.pdf")},
    content = function(file) {
      pdf(file, width=8.5, height=5)
      myplots2()
      dev.off()      
    })


})














