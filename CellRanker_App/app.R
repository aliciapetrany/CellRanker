#CellRankers R shiny application

library(shiny)
library(plotly)
library(Seurat)
library(stringi)

example.cords <- read.csv("example_tsne_coords.csv") #reads example data
outputdf <- example.cords[,c(1,4)]

# Define UI for application 
ui <- fluidPage(

    #set title
    titlePanel("CellRanker"),
    
      sidebarPanel(
        h5("Upload MTX files"), #section header
        
        fileInput("file1", "Barcodes"), # Intput file 1
        
        fileInput("file2", "Features"), # Input file 2
        
        fileInput("file3", "Matrix"), # Input file 3
       
        tags$hr(), #horizontal line
        
        h5("Or, use the example dataset"), #secion header
        
        actionButton("example", "Example"), #button to use example data
       
        tags$hr(),
          
        downloadButton("downloadData", "Download Results") #button for downloads
        
        ),
    
    mainPanel(plotlyOutput("graph")) #interactive tsne plot
)
    
    
  
# Define server logic 
server <- function(input, output) {
  options(shiny.maxRequestSize=30*1024^2) #changes max upload size
  
  #Observe block to regulate CellRanker Workflow
    observe({
      #if all files are uploaded, CellRanker begins running
        req(input$file1)
        req(input$file2)
        req(input$file3)
        showNotification("CellRanker Is Running. Take A Seat And Get Some Coffee")
        
      #Runs CellRanker workflow         
        mtx <- ReadMtx(mtx = input$file3$datapath, cells = input$file1$datapath, features = input$file2$datapath )
        User_SeuratObj <- CreateSeuratObject(mtx)
        User_SeuratObj <- CellRanker(User_SeuratObj)
        tsne_coords <- data.frame(User_SeuratObj@reductions$tsne@cell.embeddings)
        tsne_coords$labels <- User_SeuratObj@meta.data$cell_type
    
    #stores labels in a dataframe with their respective barcode for user to download
      outputdf <- data.frame(barcode = colnames(User_SeuratObj@assays$RNA@data),
               celltype = User_SeuratObj@meta.data$cell_type) 
      
     #plots the output
       output$graph <- renderPlotly({
         plot_ly(tsne_coords, 
                 x = ~tSNE_1, 
                 y = ~tSNE_2, 
                 type = 'scatter', 
                 mode = 'markers',
                 split = ~labels,
                 marker = list(size = 2))
          })
    })
    
    #Event to manage visualizing sample data in response to button press
    observeEvent(input$example, {
        output$graph <- renderPlotly({
            plot_ly(example.cords, 
                   x = ~tSNE_1, 
                   y = ~tSNE_2, 
                   type = 'scatter', 
                   mode = 'markers',
                   split = ~labels,
                   marker = list(size = 2))
         })
    })
        
    #Event to download either example or user results.
    #Output file is a csv containing barcods and label
    output$downloadData <- downloadHandler(
         filename = function() {
           paste("input$dataset", ".csv", sep = "")
         },
         content = function(file) {
           write.csv(outputdf, file)
         }
    )
  
}
    

# Run the app
shinyApp(ui = ui, server = server)




