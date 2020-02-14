
Packages <- c("dplyr","tidyr","ggplot2","ggfortify","factoextra","pheatmap","gplots","mixOmics","qgraph","corrplot","data.table","gridExtra","sna","intergraph","igraph","Matrix","ggrepel","knitr","tidyr")
invisible(lapply(Packages, library, character.only = TRUE))

ExploratoryAnalysis = function(data=NULL)
{

  
  #### -------------------> Correlation
  col3 <- colorRampPalette(c("steelblue4","white","firebrick")) 
  corrplot.mixed(cor(data[,-1], use = "complete.obs"), number.cex = .7, lower.col = 'black', upper.col = col3(100))
  
  p_corr <- recordPlot()
  
  
  #### -------------------> PCA
  pca.X <- mixOmics::pca(t(data[,-1]),center=T,scale=F)
  
  loadings_PC1 <- data.frame(mixOmics::selectVar(pca.X, comp=1)$value)
  row.names(loadings_PC1) <- data$Knockout[order(loadings_PC1)]
  loadings_PC2 <- data.frame(mixOmics::selectVar(pca.X, comp=2)$value)
  row.names(loadings_PC2) <- data$Knockout[order(loadings_PC2)]
  PCA_loadings <- data.frame(cbind(loadings_PC1,loadings_PC2))
  names(PCA_loadings) <- c('PC1','PC2')
  
  # Individual factor map
  dtp <- data.frame('Ion' = row.names(data.frame(pca.X$variates)), pca.X$variates)
  names(dtp) <- c("Ion","PC1","PC2")
  
  pca_plot <- ggplot(data = dtp, aes(x = PC1, y = PC2)) + 
    geom_point(color='steelblue', size=3, alpha=0.4) +
    geom_text_repel(aes(label = row.names(data.frame(pca.X$variates))), size=4) +
    theme_bw() +
    theme(legend.position = "none") +
    xlab(paste("PC1: ",round(pca.X$explained_variance[1]*100,2),"% expl. variance", sep='')) +
    ylab(paste("PC2: ",round(pca.X$explained_variance[2]*100,2),"% expl. variance", sep='')) +
    labs(title = "PCA")
  
  
  #### -------------------> HEATMAP
  pheatmap(data[,-1], 
           show_rownames=F, cluster_cols=T, cluster_rows=T, 
           legend = T,
           fontsize = 15,
           clustering_method = 'ward.D',
           scale="row")
  
  pheat <- recordPlot()
  
  
  #### -------------------> PAIRWISE CORRELATION MAP
  col <- colorRampPalette(c("skyblue4", "white", "plum4"))(20)
  
  corrI <- cor(na.omit(data[,-1]))
  heatmap(x = corrI, col = col, symm = TRUE, cexRow=1.4, cexCol=1.4, main = '')
  
  pcm <- recordPlot()
  
  #### -------------------> Regularized partial correlation network MAP
  cad <- cor_auto(data[,-1])
  qgraph(cad, graph = "glasso", layout = "spring", 
         tuning = 0.25,sampleSize = nrow(data[,-1]))
  
  Graph_lasso <- recordPlot()
  
  
  #### -------------------> Output
  res <- list()
  class(res) = "ExploratoryAnalysis"
  
  res$plot.Pearson_correlation <- p_corr
  res$plot.PCA_Individual <- pca_plot
  
  res$data.PCA_loadings <- PCA_loadings

  res$plot.heatmap <- pheat
  res$plot.pairwise_correlation_map <- pcm
  res$plot.regularized_partial_correlation_network <- Graph_lasso 
  
  return(res)
}


