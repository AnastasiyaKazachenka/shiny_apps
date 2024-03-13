library(shiny)
library(tidyverse)
library(survival)
library(survminer)
library(reshape2)
library(stats)
library(ggpubr)
library(regioneR)
library(karyoploteR)
library(CopyNumberPlots)

survivalminer <- function(surx, titx) {
  colnames(surx) <- c("SUBJID","Time","Event","Marker")
  event.code <- 0
  
  if (length(unique(sort(surx$Marker)))==2) {
    
    surx$Marker<- as.factor(surx$Marker)
    
    os <- with(surx, Surv(Time, Event == event.code))
    fit <- surv_fit(Surv(Time, Event == event.code) ~ Marker, data=surx, conf.type = "log-log")
    
    # hr w/ 95% CIs:
    m1 <- with(surx, coxph(os ~ Marker)) 
    m1.sum <- summary(m1)
    
    pv <- round(m1.sum$sctest["pvalue"],3)
    if(pv ==0){pval <- "p < 0.001"} else {pval <- paste0("p = ",pv)}
    
    hr <- round(m1.sum$conf.int[1],2)
    
    g <- ggsurvplot(fit, conf.int=FALSE, pval=FALSE, size = 1, palette = "jco",
                    risk.table=T,risk.table.y.text.col = T, fontsize=6,
                    risk.table.height = 0.4, surv.plot.height = 0.7,
                    ncensor.plot = FALSE, 
                    tables.theme = theme_survminer(font.main = 16, font.x = 16, font.tickslab = 16),
                    ylab = "Survival probability", 
                    break.time.by = 6,
                    font.legend = c(16,"plain", "black"),
                    font.main = c(16, "plain", "black"),
                    font.x = c(16, "plain", "black"), font.y = c(16, "plain", "black"),
                    font.tickslab = c(16, "plain", "black"),
                    legend="none",
                    title = titx,
                    xlab = "Time in months")
    
    g$plot <- g$plot +
      ggplot2::annotate(
        "text",
        x = Inf, y = Inf,
        vjust = 1, hjust = 1,
        label = paste0("HR = ",hr,"\n(95% CI ",round(m1.sum$conf.int[3],2),"-",round(m1.sum$conf.int[4],2),")","\n p = ",round(m1.sum$sctest["pvalue"],3)),
        size = 6
      ) } else {
        surx$Marker<- as.factor(surx$Marker)
        
        os <- with(surx, Surv(Time, Event == event.code))
        fit <- surv_fit(Surv(Time, Event == event.code) ~ Marker, data=surx, conf.type = "log-log")
        
        # hr w/ 95% CIs:
        m1 <- with(surx, coxph(os ~ Marker)) 
        m1.sum <- summary(m1)
        
        pv <- round(m1.sum$sctest["pvalue"],3)
        if(pv ==0){pval <- "p < 0.001"} else {pval <- paste0("p = ",pv)}
        
        hr <- round(m1.sum$conf.int[1],2)
        
        g <- ggsurvplot(fit, conf.int=FALSE, pval=FALSE, size = 1, palette = "jco",
                        risk.table=T,risk.table.y.text.col = T, fontsize=6,
                        risk.table.height = 0.4, surv.plot.height = 0.7,
                        ncensor.plot = FALSE, 
                        tables.theme = theme_survminer(font.main = 16, font.x = 16, font.tickslab = 16),
                        ylab = "Survival probability", 
                        break.time.by = 6,
                        font.legend = c(16,"plain", "black"),
                        font.main = c(16, "plain", "black"),
                        font.x = c(16, "plain", "black"), font.y = c(16, "plain", "black"),
                        font.tickslab = c(16, "plain", "black"),
                        legend="none",
                        title = titx,
                        xlab = "Time in months")
        
        g$plot <- g$plot +
          ggplot2::annotate(
            "text",
            x = Inf, y = Inf,
            vjust = 1, hjust = 1,
            label = paste0("p = ",round(m1.sum$sctest["pvalue"],3)),
            size = 6)
      }
  print(g)  
}




options(shiny.useragg = FALSE)