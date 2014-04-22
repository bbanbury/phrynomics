library(shiny)
shinyUI(pageWithSidebar(
  
  headerPanel(
    h1("Phrynomics")),

  sidebarPanel(
    a("See phrynomics GitHub for full code base", href="https://github.com/bbanbury/phrynomics"),
    br(),
    br(),
    h4("Original Dastset"),    
    textInput("datasetName", "Enter SNP Dataset Name:", value="SNPdata"), 
    fileInput("SNPdataset", "Choose File To Upload (.phy or .nex):", accept=c("text/csv", "text/comma-separated-values", "text/plain", ".txt", ".nex")),
#    a("See example input", href="https://raw.githubusercontent.com/bbanbury/phrynomics/master/shinyPhrynomics/example.txt"),
    numericInput("obs", "Number Of Taxa to Preview:", 5),
    checkboxInput("rmInvSites", "Remove Invariant Sites", value=TRUE),
    checkboxInput("rmNonBin", "Remove Non-Binary Sites", value=TRUE),
    checkboxInput("takeRandom", "Take a single SNP from each locus (original data must be separated by a space)", value=FALSE),
    br(),
    h4("Transformations"),
    selectInput("transformChar", "Character for Missing Data:", choices=c("-", "?", "N"), selected="?"),
    h5("Transform to SNAPP"),
    h6(em("This will translate only binary sites. The most common allele will be translated to a 0, the least common to a 2, and the heterozygote to a 1.")),
    checkboxInput("transSNAPP", "Translate", value=FALSE),
    br(),
    h5("Transform to MrBayes"),
    h6(em("A will always return a 1, T=2, G=3, C=4, and ambiguity codes get parentheses with possible bases")),
    checkboxInput("transMrBayes", "Translate", value=FALSE),
    br(),
    img(src="http://upload.wikimedia.org/wikipedia/commons/4/46/Horned_lizard_032507_kdh.jpg", width=350)),

  mainPanel(
    h4("Preview Original Dataset"),
    tableOutput("inputPreview"),
    h4("Preview Results"),
    tableOutput("resultsPreview"),
    verbatimTextOutput("communicationWindow"),
    downloadButton("downloadPhy", label="Download Full Results to Phy"),
    downloadButton("downloadNex", label="Download Full Results to Nexus"))
))
