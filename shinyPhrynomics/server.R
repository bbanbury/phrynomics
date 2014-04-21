library(shiny)
library(devtools)
source("~/phrynomics/trunk/phrynomicsFunctions.R")


##Communicates with UI
shinyServer(function(input, output) {
  contents <- reactive({
    inFile <- input$SNPdataset
    if (is.null(inFile))
      return(NULL)
    results <- read.table(inFile$datapath, row.names=1, skip=1, stringsAsFactors=FALSE)
    if(input$rmInvSites)
      results <- RemoveInvariantSites(results)
    if(input$rmNonBin)
      results <- RemoveNonBinary(results)
    if(input$takeRandom)
      results <- TakeSingleSNPfromEachLocus(results)[[1]]
    if(input$transSNAPP)
      results <- TranslateBases(results, ordered=TRUE)
    if(input$transMrBayes)
      results <- TranslateBases(results, ordered=FALSE)
    return(results)
  })  
    
  output$inputPreview <- renderTable({  
    inFile <- input$SNPdataset
    if (is.null(inFile))
      return(NULL)
    dataTable <- read.table(inFile$datapath, row.names=1, skip=1)[1:input$obs,]
  }, include.colnames=FALSE)
  
  output$resultsPreview <- renderTable({   
    outPrev <- head(contents(), n=input$obs)
    if(!is.null(outPrev))
      outPrev <- cSNP(SplitSNP(outPrev))
  }, include.colnames=FALSE)
  
  output$downloadData <- downloadHandler(
    filename = function() { paste(input$datasetName, ".txt", sep="") },
    content = function(file) {
      fileContent <- paste(dim(contents())[1], sum(nchar(contents()[1,])))
      write(fileContent, file)
      write.table(contents(), file, quote=FALSE, append=TRUE, col.names=FALSE)
    }
  )

  
#  myInput <- reactive({read.table("input$SNPdataset")})
#  newData <- reactive({IsVariable(input$data)})
#  output$print_newData <- renderPrint(newData())
#  output$print_myInput <- renderPrint(myInput())
})

