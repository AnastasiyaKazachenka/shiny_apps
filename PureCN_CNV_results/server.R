load("tumor_202_DNAseqments.RData")
load("tumor_102_DNAseqments.RData")
load("PureCN_results.RData")
M_202 <- c("1401001","8704001","8720001","8704002","1401003","8709001","8720002","8705003","8720003","8713001","8702001","8707001","1401005","3703001","8714003","8703001","8705005","8709002","1902001","8707002","8714004","8201002","1401006","8710001","1201001","8708002","8750002","8710004","8201003","6101001","8706002","1902002","6101002","1901004","8706001","6101004","6101005","6101006","8201004","3501002","8705006","8705008","8710005","8201005","6101007","1204001","1201002","8714005","8714007","8705009","8707003","8707004","6101008","6101009","3504006","3504007","8705010","1401010","6204001","6101010","8715004","8718003","8719004","8707007","3506004","8702003","1202002","8201007","5901001","5901002","3503002","1201004","6101014","1401012","8707005","8708005","8714008","8201009","1201005","6101016","3701001","3504010","6401001","8704004","8705013","8705011","1901010","1901005","1901011","8718004","8702004","8719005","3701012","3503003","3501008","1401011","3701007","3702002","6401006","6401005","6207001","8713006","1201006","8720007","6502001","6501001","6502002","8702006","6207005","3701020","1202005","1201007","3503004","3503006","3503008","8704005","8718007","8718006","8701001","1401015","3501010","3504014","6101021","5002001","3504015","6101025","6102005","6401018","6101023","6501002","6401015","1401016","8710009","8704006","8723001","8723002","8718010","5901003","3501011","3501012","3506011","8201011","3507002","3507004","6102007","6204003","3501013","6502006","6101029","6401011","1902004","8720009","8719007","8710010","8708007","8702009","3506010","8201010","8201012","6101032","6101033","6302007","6102008","6101034","3504017","6302005","8707010","8724002","1401017","1201009","8201013","6302001","6401023","6401022","1901015","8704010","1901016","1901017","8714010","3702007","3702010","6401020","6401028")
P_202 <- c("1202001","8201001","8713005","3503001","3501005","3703003","3501006","8714001","1401009","3703004","3506002","3504009","3504008","6401002","3701013","3701016","3501009","6207002","3701015","6207004","3701018","6401010","6401007","6101017","5001001","3701019","8721005","6401014","1202008","3508001","6302002","6302008","6302006","6302003","3703005","6401017")

tebe_202 <- c("1201001","1201002","1201003","1201005","1201006","1201009","1202002","1202004","1202008","1401001","1401002","1401004","1401006","1401007","1401010","1401011","1401012","1401013","1401015","1401017","1901001","1901003","1901004","1901005","1901010","1901011","1901013","1901014","1901015","1901016","1901017","1901018","1902001","1902004","3501002","3501006","3501008","3501009","3501012","3503001","3503004","3503006","3503007","3503008","3504001","3504002","3504003","3504004","3504005","3504007","3504008","3504013","3504014","3504015","3504016","3504018","3505001","3505002","3505003","3505004","3506002","3506004","3506006","3506007","3506009","3506010","3506012","3506014","3507004","3508002","3701001","3701002","3701003","3701007","3701009","3701010","3701012","3701014","3701015","3701016","3701017","3701018","3701019","3701020","3701022","3702001","3702002","3702005","3702007","3702009","3702010","3702011","3702012","3703003","3703004","3703005","5001001","5001003","5001004","5002001","5901001","5901002","5901005","6101002","6101005","6101006","6101007","6101010","6101011","6101012","6101013","6101014","6101015","6101016","6101018","6101021","6101022","6101023","6101026","6101028","6101034","6102002","6102005","6102007","6102008","6102010","6204001","6207001","6207003","6207004","6207005","6207006","6207010","6207011","6208001","6208003","6208004","6302001","6302002","6302003","6302004","6302005","6302006","6302007","6302008","6401001","6401002","6401004","6401007","6401009","6401010","6401011","6401014","6401017","6401022","6401023","6401024","6401026","6401027","6401028","6501001","6501002","6501003","6502002","6502004","6502005","6502006","8201001","8201002","8201004","8201005","8201006","8201011","8201013","8701001","8701002","8702001","8702003","8702004","8702005","8702009","8703001","8703003","8704002","8704005","8704008","8704009","8705004","8705005","8705006","8705008","8705009","8705011","8705014","8705015","8706001","8707002","8707003","8707006","8707008","8708001","8708003","8708005","8708006","8709001","8709002","8709003","8710004","8710005","8710009","8711001","8713001","8713004","8713005","8713007","8714001","8714003","8714004","8714005","8714008","8714010","8715002","8715004","8718002","8718003","8718004","8718006","8718007","8718009","8718010","8719004","8719006","8719007","8720002","8720004","8720005","8721007","8722001","8722003","8722004","8722005","8723001","8723002","8750001","8750003")
inv_202 <- c("1201004","1201007","1202001","1202005","1202006","1204001","1401003","1401009","1901002","1902002","3501011","3501013","3503002","3503003","3504006","3504010","3504017","3506001","3506003","3506011","3506013","3507001","3507002","3507003","3508001","3701005","3701011","3701013","3702003","3702006","3703001","5901003","5901004","6101001","6101004","6101008","6101009","6101017","6101025","6101027","6101029","6101032","6101033","6102003","6102004","6207002","6207008","6208002","6209001","6401005","6401006","6401015","6401018","6401019","6502001","6504001","8201003","8201009","8201012","8702002","8702007","8703002","8704001","8704004","8704006","8704010","8705003","8705010","8705013","8707001","8707004","8707005","8707007","8707009","8707010","8708002","8710001","8710008","8713006","8714002","8714006","8714007","8718005","8719001","8719003","8719005","8720001","8720003","8720007","8720008","8720009","8721005","8724001","8724002","8750002")

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