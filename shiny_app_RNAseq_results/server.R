mat <- readRDS("Esra_vsd_matrix_for_PCA.rds")
annot <- read.csv("esra_annotation.csv",header=T)
genes_annot <- readRDS("Genes_Ensemble_to_ID.rds")
TPMs <- readRDS("TPM_counts_with_gene_info.rds")
row.names(annot) <- annot$SampleID
PWs_1 <- readRDS("DE_pathways_KEGG_REACTOME_BIOCARTA_TPMs_Poisson.rds")
geneset_1 <- readRDS("KEGG_BIOCARTA_REACTOME_ensemblegenes.rds")
load("DE_genes_all_cell_lines.RData")
PWs_2 <- readRDS("DE_pathways_GO_TPMs_Poisson.rds")
geneset_2 <- readRDS("GO_ensemblegenes.rds")

function(input,output,session) {
  
### PAC page
  gene_list <- reactive({
    if (!(input$gene=="")) {
      genes <- strsplit(input$gene,",")[[1]]
      genes2 <- genes_annot$genes[genes_annot$new_genes %in% genes]
      return(genes2)
    }
  })
  
  my_data <- reactive(
  if (input$gene=="") {
    pcaData2 <- prcomp(mat, rank.=3)
    components <- pcaData2[["x"]]
    components <- data.frame(components)
    components$PC2 <- -components$PC2
    components$PC3 <- -components$PC3
    components <- merge(components,annot, by=0)
    row.names(components) <- components$Row.names
    components$Row.names <- NULL
    components <- components %>%
      filter(Cell_line %in% input$cell_line,Treatment %in% input$treatment)
    return(components)
  } else {
    mat2 <- mat[,gene_list()]
    pcaData2 <- prcomp(mat2, rank.=3)
    components <- pcaData2[["x"]]
    components <- data.frame(components)
    components$PC2 <- -components$PC2
    components$PC3 <- -components$PC3
    components <- merge(components,annot, by=0)
    row.names(components) <- components$Row.names
    components$Row.names <- NULL
    components <- components %>%
      filter(Cell_line %in% input$cell_line,Treatment %in% input$treatment)
    return(components)
  }
)

my_groups <- reactive(my_data()[,input$color_category])
my_colors <-reactive(get_colors(my_groups()))
legend_colors <- reactive(as.data.frame(my_colors(),row.names = names(my_colors())) %>% distinct())

  output$PCA_plot <- renderRglwidget({
    try(close3d(),silent=TRUE)
    spheres3d(my_data()$PC1,my_data()$PC2,my_data()$PC3,r=5,color = my_colors())
#    rgl_add_axes(my_data()$PC1,my_data()$PC2,my_data()$PC3,show.bbox = FALSE)
    axes3d(edges=c("x--","y--","z--"))
    title3d(xlab="PC1",ylab="PC2",zlab="PC3")
    rglwidget()
  })
  
  output$Legend <- renderPlot({
    plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
    legend("topleft", legend =row.names(legend_colors()), pch=16, pt.cex=3, cex=1.5, bty='n',
           col =legend_colors()[,1])
    mtext(input$color_category, at=0.2, cex=2)
  })
  
## Heatmap page

  output$Heatmap <- renderPlot({
    if (input$gene=="") {
      my_plot <- ggplot() + annotate("text", x = 10,  y = 10, size = 6, label = "O       o O       o O       o
| O   o | | O   o | | O   o |
| | O | | | | O | | | | O | |
| o   O | | o   O | | o   O |
o       O o       O o       O

I need at least 2 genes
") + theme_void()
    } else {
      annot2 <- annot %>%
        filter(Cell_line %in% input$cell_line,Treatment %in% input$treatment) %>%
        arrange(Cell_line,Treatment)
      TPMs2 <- TPMs[,annot2$SampleID]
      TPMs2 <- merge(genes_annot,TPMs2,by=0)
      TPMs2 <- TPMs2[rowSums(TPMs2[,annot2$SampleID]) > 0,]
      if (input$log == "no") {
        TPMs2 <- TPMs2 %>%
          filter(genes %in% gene_list()) %>%
          na.omit() %>%
          mutate(new_genes=make.unique(new_genes)) %>%
          set_rownames(.$new_genes) %>%
          dplyr::select(annot2$SampleID)
      } else {
        TPMs2 <- TPMs2 %>%
        filter(genes %in% gene_list()) %>%
        na.omit() %>%
        mutate(new_genes=make.unique(new_genes)) %>%
        set_rownames(.$new_genes) %>%
        dplyr::select(annot2$SampleID) %>%
        mutate_if(is.numeric,log1p)
      }
      my_plot <- pheatmap(TPMs2,scale=input$scaled,annotation_col = annot2[, c("Treatment","Cell_line")],cluster_cols=as.logical(input$cluster_cols))
    }
    my_plot
    })
  
##Boxplot page
  output$Boxplot1 <- renderPlot({
    if ((input$gene=="")) {
      my_plot2 <- ggplot() + annotate("text", x = 10,  y = 10, size = 6, label = "O       o O       o O       o
| O   o | | O   o | | O   o |
| | O | | | | O | | | | O | |
| o   O | | o   O | | o   O |
o       O o       O o       O
                                      
Too little ... I need 1 gene") + theme_void()
    } else if (length(gene_list()) >1 ) {
      my_plot2 <- ggplot() + annotate("text", x = 10,  y = 10, size = 6, label = "O       o O       o O       o
| O   o | | O   o | | O   o |
| | O | | | | O | | | | O | |
| o   O | | o   O | | o   O |
o       O o       O o       O
                                      
Too many ... I need 1 gene") + theme_void() } else {
      my_plot2 <- TPMs %>%
        pivot_longer(!c(genes,new_genes),names_to="SampleID",values_to = "TPMs") %>%
        inner_join(annot) %>%
        filter(genes %in% gene_list()) %>%
        filter(Cell_line %in% input$cell_line) %>%
        filter(Treatment %in% input$treatment) %>%
        ggplot(aes(x=factor(Cell_line,level=c("A375","MP41","OV56","SNU685","NCIH441","SW620","MEL624","NCIH1755")),y=TPMs,color=Treatment)) +
        geom_boxplot(position=position_dodge(0.9),outlier.shape = NA) +
        geom_dotplot(binaxis = 'y',stackdir='center',position=position_dodge(0.9),dotsize=0.5) +
        ggtitle(input$gene) +
        geom_hline(yintercept = 5, color="red") +
        theme_bw()
    }
    my_plot2
  })
  
## comparison page
  
  output$Boxplot2 <- renderPlot({
    if ((input$gene=="" | length(gene_list())==1)) {
      my_plot3 <- ggplot() + annotate("text", x = 10,  y = 10, size = 6, label = "O       o O       o O       o
| O   o | | O   o | | O   o |
| | O | | | | O | | | | O | |
| o   O | | o   O | | o   O |
o       O o       O o       O
                                      
Too little ... I need 2 genes") + theme_void()
    } else if (length(gene_list()) >2 ) {
      my_plot3 <- ggplot() + annotate("text", x = 10,  y = 10, size = 6, label = "O       o O       o O       o
| O   o | | O   o | | O   o |
| | O | | | | O | | | | O | |
| o   O | | o   O | | o   O |
o       O o       O o       O
                                      
Too many ... I need 2 genes") + theme_void() } else {
  my_tpms <- TPMs %>%
    filter(genes %in% gene_list()) %>%
    pivot_longer(!c(genes,new_genes),names_to="SampleID",values_to = "TPMs") %>%
    dplyr::select(-genes) %>%
    pivot_wider(names_from = "new_genes",values_from = "TPMs") %>%
    inner_join(annot) %>%
    filter(Cell_line %in% input$cell_line) %>%
    filter(Treatment %in% input$treatment)
  genes_to_plot <- colnames(my_tpms[,c(2,3)])
  my_plot3 <- my_tpms %>%
    f1(genes_to_plot)
}
    my_plot3
  })
  
  ## Pathways page
  
  PW_DE <- reactive({
    if (input$pws_collection == "KEGG,REACTOME,BIOCARTA") {
      minset_subset<- subEset(PWs_1,subset=list(Treatment=c("Control supernatant","ImmTAC supernatant"),Cell_line=input$pws))
      mod <- model.matrix(~ factor(minset_subset$Treatment))
      colnames(mod) <- c("all","controlvsImmTAC")
      fit <- lmFit(minset_subset,mod)
      fit <- eBayes(fit)
    } else {
      minset_subset<- subEset(PWs_2,subset=list(Treatment=c("Control supernatant","ImmTAC supernatant"),Cell_line=input$pws))
      mod <- model.matrix(~ factor(minset_subset$Treatment))
      colnames(mod) <- c("all","controlvsImmTAC")
      fit <- lmFit(minset_subset,mod)
      fit <- eBayes(fit)
    }
    return(fit)
  })
  
  output$Volcano <- renderPlot({
    tt <- topTable(PW_DE(), coef=2, n=Inf)
    DEpwys <- rownames(tt)[tt$adj.P.Val <= 0.05]
    plot(tt$logFC, -log10(tt$P.Value), pch=".", cex=4, col=grey(0.75),
         main=input$pws, xlab="GSVA enrichment score difference", ylab=expression(-log[10]~~Raw~P-value))
    abline(h=-log10(max(tt$P.Value[tt$adj.P.Val <= 0.05])), col=grey(0.5), lwd=1, lty=2)
    points(tt$logFC[match(DEpwys, rownames(tt))],
           -log10(tt$P.Value[match(DEpwys, rownames(tt))]), pch=".", cex=5, col="darkred")
    text(max(tt$logFC)*0.85, -log10(max(tt$P.Value[tt$adj.P.Val <= 0.05])), "5% FDR", pos=3)
  })
  
  output$PW_heatmap <- renderPlot({
    tt <- topTable(PW_DE(), coef=2, n=Inf)
    DEpwys <- rownames(tt)[tt$adj.P.Val <= 0.05]
    annot3 <- annot %>%
      filter(Cell_line %in% input$cell_line,Treatment %in% input$treatment)
    if (input$pws_collection == "KEGG,REACTOME,BIOCARTA") {
      PWs_red <- exprs(PWs_1[DEpwys,])
    } else {
      PWs_red <- exprs(PWs_2[DEpwys,])
    }
    my_plot <- PWs_red[,annot3$SampleID] %>%
      pheatmap(scale=input$scaled,annotation_col = annot3[, c("Treatment","Cell_line")],show_rownames = F)
  })
  
  PW_table_v <- reactive({
    tt <- topTable(PW_DE(), coef=2, n=Inf)
    tt <- tt[tt$adj.P.Val <= 0.05,]
    if (input$pws_collection == "KEGG,REACTOME,BIOCARTA") {
      tt <- merge(geneset_1,tt,by=0) %>%
        dplyr::select(-logFC,-AveExpr,-t,-P.Value,-B) %>%
        arrange(adj.P.Val)
      row.names(tt) <- tt$Row.names
      tt$Row.names<-NULL
    } else {
      tt <- merge(geneset_2,tt,by=0) %>%
        dplyr::select(-logFC,-AveExpr,-t,-P.Value,-B) %>%
        arrange(adj.P.Val)
      row.names(tt) <- tt$Row.names
      tt$Row.names<-NULL
    }
    tt})
  output$PW_table <- renderDataTable({
    PW_table_v()
  })
  

### DE genes page
  
  deGenes <- reactive({
    deGenes_c <- get(input$pws) %>%
      filter(padj < 0.05) %>%
      filter(abs(log2FoldChange) >= 1) %>%
      select(genes)
    return(deGenes_c$genes)
  })
  geneUniverse <- reactive({
    geneUniverse_c <- get(input$pws) %>%
      dplyr::select(genes)
    return(geneUniverse_c$genes)
  })
  
  output$GeneDE_heatmap <- renderPlot({
    annot2 <- annot %>%
      filter(Cell_line %in% input$cell_line,Treatment %in% input$treatment) %>%
      arrange(Cell_line,Treatment)
    TPMs2 <- TPMs[,annot2$SampleID]
    TPMs2 <- merge(genes_annot,TPMs2,by=0)
    TPMs2 <- TPMs2[rowSums(TPMs2[,annot2$SampleID]) > 0,]
    if (input$log == "no") {
      TPMs2 <- TPMs2 %>%
        filter(genes %in% deGenes()) %>%
        na.omit() %>%
        mutate(new_genes=make.unique(new_genes)) %>%
        set_rownames(.$new_genes) %>%
        dplyr::select(annot2$SampleID)
    } else {
      TPMs2 <- TPMs2 %>%
        filter(genes %in% deGenes()) %>%
        na.omit() %>%
        mutate(new_genes=make.unique(new_genes)) %>%
        set_rownames(.$new_genes) %>%
        dplyr::select(annot2$SampleID) %>%
        mutate_if(is.numeric,log1p)
    }
    my_plot <- pheatmap(TPMs2,scale=input$scaled,annotation_col = annot2[, c("Treatment","Cell_line")],cluster_cols=as.logical(input$cluster_cols),show_rownames=F,main=paste(nrow(TPMs2)," DE genes for ",input$pws, "Control SN vs ImmTAC SN"))
    my_plot
  })
  
  
  
  output$GO_plot <- renderPlot({
    gene_num <- length(deGenes())
    deGenes_c2 <- unlist(mget(deGenes(),envir=org.Hs.egENSEMBL2EG,ifnotfound=NA))
    geneUniverse_c2 <- unlist(mget(geneUniverse(),envir=org.Hs.egENSEMBL2EG,ifnotfound=NA))
    ans.go <- enrichGO(gene=deGenes_c2, ont="BP",OrgDb = "org.Hs.eg.db",universe=geneUniverse_c2,readable = TRUE, pvalueCutoff=0.05)
    go_plot <- dotplot(ans.go, showCategory=20) + ggtitle(paste("Enrichment of GO terms among ", gene_num," DE genes for ",input$pws))
    go_plot
  })
  
  output$GeneDE_table <- renderDataTable({
    DE_table <- get(input$pws) %>%
      filter(padj < 0.05) %>%
      filter(abs(log2FoldChange) >= 1)
    DE_table
  })
  

##Download tables
  my_tpms <- reactive({
    my_tpms <- TPMs %>%
    filter(genes %in% gene_list()) %>%
    pivot_longer(!c(genes,new_genes),names_to="SampleID",values_to = "TPMs") %>%
    dplyr::select(-genes) %>%
    pivot_wider(names_from = "new_genes",values_from = "TPMs") %>%
    inner_join(annot) %>%
    filter(Cell_line %in% input$cell_line) %>%
    filter(Treatment %in% input$treatment)
    
    return(my_tpms)})
                          
  output$downloaddata <- downloadHandler(filename = "TPMs_table.csv",
                                         content = function(file) {
                                           write.csv(my_tpms(),file,row.names = TRUE)
                                         })
  
  output$downloaddata2 <- downloadHandler(filename = "DE_pathways_table.csv",
                                          content = function(file) {
                                            write.csv(PW_table_v(),file,row.names = TRUE)
                                          })
}

