load("tumor_202_DNAseqments.RData")
load("tumor_102_DNAseqments.RData")
load("PureCN_results.RData")

M_202 <- c(1:10) #should be IDs of metastatic patients
P_202 <- c(11:20) #should be IDs of primary patients

tebe_202 <- c(1:10) #should be ids of patients with tebe treatment
inv_202 <- c(11:20) #should be ids of investigator choice treatment

genes <- read.csv("./genes_coords.csv")
annot <- read.csv("anot.csv",colClasses = c("SUBJID"="character"))
hg38 <- filterChromosomes(getGenome("hg38"), organism = "hg", chr.type = "canonical")

function(input,output,session) {
  updateSelectizeInput(session,"gene",choices=sort(unique(genes$gene.symbol)),server=TRUE,selected="BAP1")
  updateSelectizeInput(session,"chr",choices=sort(unique(genes$chr)),server=TRUE,selected="chr3")
  
  new_keys <- reactive({
    new_keys <- keys_202
    if (input$MetPr202 == "metastatic") {
      new_keys <- new_keys[new_keys %in% M_202]
    }
    if (input$MetPr202 == "primary") {
      new_keys <- new_keys[new_keys %in% P_202]
    }
    if (input$TebeInv202 == "tebe") {
      new_keys <- new_keys[new_keys %in% tebe_202]
    }
    if (input$TebeInv202 == "inverstigator choice") {
      new_keys <- new_keys[new_keys %in% inv_202]
    }
    new_keys
  })
  
  output$CN_102 <- renderPlot({
    kp <- plotKaryotype("hg38", plot.type = 4,chromosome=input$chr, cex=2)
    kpAddMainTitle(kp,main="102 dataset")
    kpAddLabels(kp,labels="log2",cex=2,srt=90,pos=3,label.margin=0.03)
    for (x in 1:length(keys_102)) {
      plotCopyNumberCallsAsLines(kp, cn.calls = get(keys_102[x]),cn.column="seg.mean",labels="",col="black",style="segments",ymax=3,ymin=-3)
    }
    kpRect(kp,chr=genes$chr[genes$gene.symbol==input$gene],x0=genes$start[genes$gene.symbol==input$gene],x1=genes$end[genes$gene.symbol==input$gene],y0=0,y1=1,data.panel="all",col="red",border="red")
  })
  
  output$CN_202 <- renderPlot({
    kp <- plotKaryotype("hg38", plot.type = 4,chromosome=input$chr, cex=2)
    kpAddMainTitle(kp,main="202 dataset")
    kpAddLabels(kp,labels="log2",cex=2,srt=90,pos=3,label.margin=0.03)
    for (x in 1:length(new_keys())) {
      if (new_keys()[x] %in% M_202) {
        plotCopyNumberCallsAsLines(kp, cn.calls = get(new_keys()[x]),cn.column="seg.mean",labels="",col="#003B7F",style="segments",ymax=3,ymin=-3)
      } else {plotCopyNumberCallsAsLines(kp, cn.calls = get(new_keys()[x]),cn.column="seg.mean",labels="",col="#FDBD13",style="segments",ymax=3,ymin=-3)}
    }
    kpRect(kp,chr=genes$chr[genes$gene.symbol==input$gene],x0=genes$start[genes$gene.symbol==input$gene],x1=genes$end[genes$gene.symbol==input$gene],y0=0,y1=1,data.panel="all",col="red",border="red")
  })
  
  output$Surv_CN_102 <- renderPlot({
    survivalminer(purecn_c_102[,c("SUBJID","AVAL","CNSR",input$gene)],paste("Based on PureCN estimates ",input$gene))
  })
  
  output$Surv_log_102 <- renderPlot({
    purecn_log <- purecn_log2_102 %>% 
      dplyr::select(SUBJID,AVAL,CNSR,input$gene)
    colnames(purecn_log) <- c("SUBJID","AVAL","CNSR","gene")
    purecn_log <- purecn_log %>% 
      mutate(gene = as.numeric(gene)) %>%
      mutate(Surv_group=case_when(
        gene > input$gain ~ "gain",
        gene <= input$gain & gene >= input$loss ~ "diploid",
        gene < input$loss ~ "loss")) %>% 
      dplyr::select(SUBJID,AVAL,CNSR,Surv_group)
    survivalminer(purecn_log,paste("Based on input thresholds ",input$gene))
  })
  
  output$log_waterfall_102 <- renderPlot({
    purecn_log <- purecn_log2_102 %>% 
      dplyr::select(SUBJID,AVAL,CNSR,input$gene)
    colnames(purecn_log) <- c("SUBJID","AVAL","CNSR","gene")
    purecn_log %>% 
      mutate(gene = as.numeric(gene)) %>%
      arrange(gene) %>%
      ggplot(aes(x=factor(SUBJID,levels=SUBJID),y=gene)) + 
      geom_bar(stat="identity") + 
      ylab(label="log2") +
      labs(x="") +
      geom_hline(yintercept = input$gain,linetype="dashed",color="red") +
      geom_hline(yintercept = input$loss,linetype="dashed",color="red") +
      theme(axis.text.x=element_blank())
  })
  
  output$Surv_CN_202 <- renderPlot({
    survivalminer(purecn_c_202[purecn_c_202$SUBJID %in% new_keys(),c("SUBJID","AVAL","CNSR",input$gene)],paste("Based on PureCN estimates ",input$gene))
  })
  
  output$Surv_log_202 <- renderPlot({
    purecn_log <- purecn_log2_202 %>% 
      dplyr::select(SUBJID,AVAL,CNSR,input$gene) %>%
      filter(SUBJID %in% new_keys())
    colnames(purecn_log) <- c("SUBJID","AVAL","CNSR","gene")
    purecn_log <- purecn_log %>% 
      mutate(gene = as.numeric(gene)) %>%
      mutate(Surv_group=case_when(
        gene > input$gain ~ "gain",
        gene <= input$gain & gene >= input$loss ~ "diploid",
        gene < input$loss ~ "loss")) %>% 
      dplyr::select(SUBJID,AVAL,CNSR,Surv_group)
    survivalminer(purecn_log,paste("Based on input thresholds ",input$gene))
  })
  
  output$log_waterfall_202 <- renderPlot({
    purecn_log <- purecn_log2_202 %>% 
      dplyr::select(SUBJID,AVAL,CNSR,input$gene) %>%
      filter(SUBJID %in% new_keys())
    colnames(purecn_log) <- c("SUBJID","AVAL","CNSR","gene")
    purecn_log %>% 
      mutate(gene = as.numeric(gene)) %>%
      arrange(gene) %>%
      inner_join(annot) %>%
      ggplot(aes(x=factor(SUBJID,levels=SUBJID),y=gene,fill=UV_type)) + 
      geom_bar(stat="identity") + 
      ylab(label="log2") +
      labs(x="") +
      geom_hline(yintercept = input$gain,linetype="dashed",color="red") +
      geom_hline(yintercept = input$loss,linetype="dashed",color="red") +
      theme(axis.text.x=element_blank())
  })
}
