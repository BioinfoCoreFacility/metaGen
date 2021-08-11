#' A cor.OTUvsGENES function
#'
#' This function allows you to calculate the correlation coeficient with p-value between two tabular data
#' @OTUs OBLIGATORY - is the dataframe with OTUs abundance
#' @genes OBLIGATORY - is the dataframe with genes abundance
#' @r is the correlation coeficient theshold
#' @p is the p-value treshold
#' @type "pearson" or "spearman" method
#' cor.OTUvsGENES(OTUs, genes)
cor.OTUvsGENES <- function(OTUs, genes, r=0.9, p = 0.05, type = "pearson"){
  require("Hmisc")
  require("reshape2")
  require("qvalue")
  OTUs <- OTUs[,sort(colnames(OTUs))]
  genes <- genes[,sort(colnames(genes))]
  if (identical(colnames(genes), colnames(OTUs))) {
    OTU.genes <- rbind(OTUs, genes)
    
    corel <- rcorr(t(as.matrix(OTU.genes)), type = type)
    cor1 <- melt(corel$r)
    cor1$p <- melt(corel$P)$value
    
    qval <- qvalue(cor1$p)
    cor1$q <- qval$qvalues
    cor1$FDR <- qval$lfdr
    
    colnames(cor1) <- c("OTUs", "genes", "r", "p", "q", "FDR")
    cor1 <- cor1[cor1$OTUs %in% rownames(OTUs),]
    cor1 <- cor1[cor1$genes %in% rownames(genes),]
    
    corel <- cor1[cor1$r > r & cor1$p < p,]
  }else{
    message("columns of both tables are not identical!!!!")
  }
  
