#' Exploratory Analysis
#'
#'Exploratory tools for Ionomics data
#' @param \code{data} a processed data frame of ion's concentrations (columns) for each gene's knockout (rows) (wide format) e.g. \code{data.wide} from the \code{PreProcessing_fn}
#'
#' @return \code{plot.Pearson_correlation} plot of Pearson's correlations across ions
#' @return \code{plot.PCA_Individual} plot of ion's profiles in a lower-dimensional space from PCA
#' @return \code{plot.heatmap} plot of ion's concentrations for each gene's knockout with added clustering dendogram
#' @return \code{plot.pairwise_correlation_map} plot of pairwise correlation coefficients across ions
#' @return \code{plot.regularized_partial_correlation_network} plot of partial correlation coefficients for pairs of ions after conditioning for the remaining ions. Red for negative, green for positive partial correlations. No connections drawn when partial correlations equals 0
#' @return \code{data.PCA_loadings} a data frame of weights of each of the original gene knockouts also known as loading vectors associated to each PC
#'
#' @examples \code{
#' library(IonFlow)
#'
#' ### Run PreProcessing function
#' pre_proc <- PreProcessing(data=IonData,stdev=pre_defined_sd)
#'
#' ### Run ExploratoryAnalysis function
#' exp_anal <- ExploratoryAnalysis(data=pre_proc$data.wide)
#'
#' # plots
#' exp_anal$plot.Pearson_correlation
#' exp_anal$plot.PCA_Individual
#' exp_anal$plot.heatmap
#' exp_anal$plot.pairwise_correlation_map
#' exp_anal$plot.regularized_partial_correlation_network
#'
#' # data
#' head(exp_anal$data.PCA_loadings)
#' }
#' @export
ExploratoryAnalysis = function(data=NULL)
{
  Packages <- c("dplyr","tidyr","ggplot2","ggfortify","factoextra","pheatmap","gplots","mixOmics","qgraph","corrplot","data.table","gridExtra","sna","intergraph","igraph","Matrix","ggrepel","knitr","tidyr")
  suppressWarnings(invisible(lapply(Packages, library, character.only = TRUE)))

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

  pca_plot <- ggplot2::ggplot(data = dtp, ggplot2::aes(x = PC1, y = PC2)) +
    geom_point(color='steelblue', size=3, alpha=0.4) +
    geom_text_repel(ggplot2::aes(label = row.names(data.frame(pca.X$variates))), size=4) +
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
  suppressWarnings(qgraph(cad, graph = "glasso", layout = "spring",
         tuning = 0.25,sampleSize = nrow(data[,-1])))

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


