fluidPage (
  sidebarLayout(
    sidebarPanel(
      selectizeInput("chr","What chromosomes to plot?",choices=NULL,selected=NULL,multiple=TRUE),
      radioButtons("MetPr202","What 202 samples to show?",c("all","metastatic","primary"),selected=c("all")),
      radioButtons("TebeInv202","What 202 treatmant to show?",c("all","tebe","inverstigator choice"),selected=c("all")),
      selectizeInput("gene","Gene for PureCN results and to be shown on CN plot",choices=NULL,selected=NULL,multiple=FALSE),
      numericInput("gain","Copy number log2 ration for GAIN",value=0.2,min=0,max=2),
      numericInput("loss","Copy number log2 ration for LOSS",value=-0.25,min=-1.5,max=0)
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("CNV plot",
                 fluidRow(plotOutput("CN_102")),
                 fluidRow(plotOutput("CN_202"))),
        tabPanel("PureCN results 102", 
                 fluidRow(plotOutput("log_waterfall_102")),
                 fluidRow(column(6,plotOutput("Surv_log_102")),
                          column(6, plotOutput("Surv_CN_102")))),
        tabPanel("PureCN results 202", 
                 fluidRow(plotOutput("log_waterfall_202")),
                 fluidRow(column(6,plotOutput("Surv_log_202")),
                          column(6, plotOutput("Surv_CN_202"))))
      )
    )
  )
)
    