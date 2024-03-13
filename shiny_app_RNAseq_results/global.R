library(shiny)
library(tidyverse)
library(ggforce)
library(plotly)
library(rgl)
library(shinyRGL)
library(ggplot2)
library(pheatmap)
library(magrittr)
library(GSVA)
library(GSEABase)
library(GSVAdata)
library(org.Hs.eg.db)
library(tidyr)
library(GOexpress)
library(limma)
library(DT)
library(clusterProfiler)
library(enrichplot)


# x, y, z : numeric vectors corresponding to
#  the coordinates of points
# axis.col : axis colors
# xlab, ylab, zlab: axis labels
# show.plane : add axis planes
# show.bbox : add the bounding box decoration
# bbox.col: the bounding box colors. The first color is the
# the background color; the second color is the color of tick marks
rgl_add_axes <- function(x, y, z, axis.col = "grey",
                         xlab = "", ylab="", zlab="",
                         show.bbox = FALSE, bbox.col = c("#333377","black"))
{ 
  
  lim <- function(x){c(-max(abs(x)), max(abs(x))) * 1.1}
  # Add axes
  xlim <- lim(x); ylim <- lim(y); zlim <- lim(z)
  rgl.lines(xlim, c(0, 0), c(0, 0), color = axis.col)
  rgl.lines(c(0, 0), ylim, c(0, 0), color = axis.col)
  rgl.lines(c(0, 0), c(0, 0), zlim, color = axis.col)
  
  # Add a point at the end of each axes to specify the direction
  axes <- rbind(c(xlim[2], 0, 0), c(0, ylim[2], 0), 
                c(0, 0, zlim[2]))
  rgl.points(axes, color = axis.col, size = 3)
  
  # Add axis labels
  rgl.texts(axes, text = c(xlab, ylab, zlab), color = axis.col,
            adj = c(0.5, -0.8), size = 2)
  
  # Add bounding box decoration
  if(show.bbox){
    rgl.bbox(color=c(bbox.col[1],bbox.col[2]), alpha = 0.5, 
             emission=bbox.col[1], specular=bbox.col[1], shininess=5, 
             xlen = 3, ylen = 3, zlen = 3) 
  }
}

#' Get colors for the different levels of 
#' a factor variable
#' 
#' @param groups a factor variable containing the groups
#'  of observations
#' @param colors a vector containing the names of 
#   the default colors to be used
get_colors <- function(groups, group.col = palette()){
  groups <- as.factor(groups)
  ngrps <- length(levels(groups))
  if(ngrps > length(group.col)) 
    group.col <- rep(group.col, ngrps)
  color <- group.col[as.numeric(groups)]
  names(color) <- as.vector(groups)
  return(color)
}

f1 <- function(df, gene_list) {
  gg <- ggplot(df, 
               aes(x = .data[[gene_list[1]]], 
                   y = .data[[gene_list[2]]], 
                   shape=Treatment)) +
    geom_point(aes(color=Cell_line),size=3) +
    geom_smooth(method=lm,se=FALSE,fullrange=TRUE,aes(linetype=Treatment),color="grey") +
    theme_bw()
  return(gg)
}

options(shiny.useragg = FALSE)
options(rgl.useNULL = TRUE)