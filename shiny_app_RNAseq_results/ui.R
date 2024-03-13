fluidPage (
  sidebarLayout(
    sidebarPanel(
      textInput("gene","Gene list",value=""),
      selectInput("color_category","Color by", c("Treatment","Cell_line","Tissue","PRAME_expression","Afffected_by_supernatant")),
      radioButtons("scaled","Heatmap, scale by",c("none","row"),selected=c("row")),
      radioButtons("cluster_cols","Heatmap cluster columns",c("yes"=TRUE,"no"=FALSE),selected = c("yes"=TRUE)),
      radioButtons("log","Log2 for Heatmap", c("yes","no"),selected = c("yes")),
      checkboxGroupInput("cell_line","What cell Line to show?",choices=c("A375","MEL624","MP41","NCIH1755","NCIH441","OV56","SNU685","SW620"),selected=c("A375","MEL624","MP41","NCIH1755","NCIH441","OV56","SNU685","SW620")),
      checkboxGroupInput("treatment","What treatment to show?",choices = c("Control supernatant","ImmTAC supernatant","Medium control"),selected = c("Control supernatant","ImmTAC supernatant","Medium control")),
      selectInput("pws","Cell line for gene and pathways DE analysis", c("A375","MEL624","MP41","NCIH1755","NCIH441","OV56","SNU685","SW620")),
      selectInput("pws_collection","Gene set collection", c("KEGG,REACTOME,BIOCARTA","GO_biological_processes")),
      downloadButton("downloaddata","Download TPMs table"),
      downloadButton("downloaddata2","Download Pathways table")
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("PCA plot",
                fluidRow(column(6,rglwidgetOutput("PCA_plot")),
                         column(6,plotOutput("Legend")))),
        tabPanel("TPMs Heatmap", plotOutput("Heatmap")),
        tabPanel("TPMs Boxplots",plotOutput("Boxplot1")),
        tabPanel("TPMs Comparison",plotOutput("Boxplot2")),
        tabPanel("DE analysis of pathways",plotOutput("Volcano"),plotOutput("PW_heatmap"),dataTableOutput("PW_table")),
        tabPanel("DE analysis of genes",plotOutput("GeneDE_heatmap"), plotOutput("GO_plot"),dataTableOutput("GeneDE_table"))
    )
  )
)
)