#' Network
#'
#'Network tools for Ionomics data
#' @param \code{data} a processed data frame of ion's concentrations (columns) for each gene's knockout (rows) (wide format) e.g. \code{data.wide} from the \code{PreProcessing_fn}
#' @param \code{data_Symb} a processed data frame of symbolised (values -1,0, or 1) ion's concentration profiles (columns) for each gene's knockout (rows) (wide format) e.g. \code{data.wide_Symb} from the \code{PreProcessing_fn}
#'
#' @return \code{stats.impact_betweeness} statistics of the impact betweenees data
#' @return \code{stats.impact_betweeness_by_cluster} statistics of the impact betweenees data by cluster
#'
#' @return \code{plot.pnet plot} of the gene's network
#' @return \code{plot.impact_betweenees} plot of the impact betweenees
#'
#' @examples \code{
#' library(IonFlow)
#'
#' ### Run PreProcessing function
#' pre_proc <- PreProcessing(data=IonData,stdev=pre_defined_sd)
#'
#' ### Run GeneNetwork function
#' gene_net <- GeneNetwork(data=pre_proc$data.wide, data_Symb=pre_proc$data.wide_Symb)
#'
#' # stats
#' gene_net$stats.impact_betweeness
#' gene_net$stats.impact_betweeness_by_cluster
#'
#' # plots
#' gene_net$plot.pnet
#' gene_net$plot.impact_betweenees
#' }
#' @export
GeneNetwork = function(data=NULL,
                       data_Symb=NULL)
{
  Packages <- c("ggrepel","ggplot2","network","intergraph","dplyr","tidyr","reshape2","reshape","ggplot2","ggfortify","Matrix","ggrepel","knitr","tidyr","GGally","factoextra","grDevices","data.table","gridExtra")
  suppressWarnings(invisible(lapply(Packages, library, character.only = TRUE)))

  # Cluster of gene with same profile
  res.dist <- dist(data_Symb[,-1], method = "manhattan")
  res.hc <- hclust(d = res.dist, method = "single")
  data_Symb$cluster <- cutree(res.hc, h = 0) # distance 0

  df <- as.data.frame(table(data_Symb$cluster)); names(df) <- c('cluster', 'nGenes')
  for(i in 1:dim(df)[1]){
    tmp_wide <- data[data_Symb$cluster==df$cluster[i],]
    label <- paste("Cluster",df$cluster[i],paste("(",df$nGenes[i], " genes)", sep=''), sep=' ')
    data_Symb$Cluster[data_Symb$cluster==df$cluster[i]] <- rep(label,nrow(tmp_wide))
  }

  # Freq in each cluster
  Rk <- as.data.frame(table(data_Symb$Cluster))
  Rk <- Rk[with(Rk, order(-Rk$Freq)), ]

  # Compute gene cluster id (index)
  index = data_Symb$Cluster

  # Cluster 2 (largest cluster) contains genes with no phenotype hence not considered (input to 0)
  index[index==Rk$Var1[1]]=0
  # Smaller clusters (Freq<10) input to 0 as well (we consider only cluster with at least 10 genes)
  for(i in which(Rk$Freq<10)[1]:dim(Rk)[1]){
    index[index==Rk$Var1[i]]=0
  }
  # Apply the cluster filtering
  index=index>0

  # cluster labels with info of accumulation/decumulation of Ions (high/lower abundance)
  x=data_Symb$Cluster[index]
  ux=unique(x)
  df.symb <- data_Symb[data_Symb$Cluster %in% c(unique(x)),]
  df <- as.data.frame(table(df.symb$Cluster)); names(df) <- c('Cluster', 'nGenes')

  # Assign label
  labelsC2 <-labelsC <- list()
  name.label <- c()
  for(i in 1:length(ux)){
    sub_df.symb <- df.symb[(df.symb$Cluster==ux[i]), 2:15]#; rownames(sub_df.symb) <- df.symb[df.symb$cluster==ux[i],1]
    if (length(names(which(colSums(sub_df.symb == 1) > 0)))>0) {l1 <- paste0(names(which(colSums(sub_df.symb == 1) > 0)),"(+)")}
    if (length(names(which(colSums(sub_df.symb == -1) > 0)))>0) {l2 <- paste0(names(which(colSums(sub_df.symb == -1) > 0)),"(-)")}
    if (length(names(which(colSums(sub_df.symb == 1) > 0)))>0) {labelsC[[i]] <- l1}
    if (length(names(which(colSums(sub_df.symb == -1) > 0)))>0) {labelsC[[i]] <- l2}
    if ((length(names(which(colSums(sub_df.symb == 1) > 0)))>0) & (length(names(which(colSums(sub_df.symb == -1) > 0)))>0)) {labelsC[[i]] <- c(l1,l2)}
    labelsC2[[i]] <- do.call(paste, c(as.list(labelsC[[i]]), sep = ", "))
  }
  names(labelsC2) = ux

  ny <- paste(ux,unlist(labelsC2),sep = ": ")
  for(i in 1:length(ux)){
    sub_df.symb <- df.symb[(df.symb$Cluster==ux[i]),]
    df.symb$Label[(df.symb$Cluster==ux[i])] <- rep(ny[i],nrow(sub_df.symb))
  }

  x1=df.symb$Label
  ux1=unique(x1)
  cpy= rainbow(length(ux))
  names(cpy) = ux1

  # Compute empirical correlation matrix
  corrGenes <- cor(t(as.matrix(data[,-1])), method = "pearson", use = "pairwise.complete.obs")

  # Subset correlation matrix based on the cluster filtering
  A = corrGenes[index,index]
  # Diagonal value (1's) put to 0 to avoid showing edges from/to the same gene
  diag(A) <- 0
  # Subset correlation matrix based on threshold=0.6
  A=(A>0.6)
  A <- ifelse(A==TRUE,1,0)
  # Generate network
  Net = network::network(A, directed = FALSE)
  #h <- network.copy(Net)
  #summary(h)

  # Network plot

  gNEt <- function(dat=NULL ){
    GGally::ggnet2(dat, color = x1, color.legend = 'Label', palette = cpy,
                   edge.alpha = 0.5, size = 2,
                   legend.size = 10, legend.position = "right")
  }


  # Impact and betweenness
  df1 <- data[data_Symb$Cluster %in% ux,] # df of 269 genes, 10 clusters
  btw <- sna::betweenness(A)
  impact <- apply(df1[,-1],1, function(a) dist(rbind(a,rep(0,ncol(df1[-1]))))) # L2 norm

  df.res <- data.frame(
    Knockout = data$Knockout[index],
    impact = round(impact,3),
    betweenness = round(btw,3),
    log.betweenness = round(log(btw+1),3),
    pos = factor(ifelse((impact < quantile(impact,.75)) & (log(btw+1) < quantile(log(btw+1),.75)),1,
                        ifelse((impact < quantile(impact,.75)) & (log(btw+1) > quantile(log(btw+1),.75)),2,
                               ifelse((impact > quantile(impact,.75)) & (log(btw+1) < quantile(log(btw+1),.75)),3,4)))),
    pos.label = factor(ifelse((impact < quantile(impact,.75)) & (log(btw+1) < quantile(log(btw+1),.75)), 'Low impact, low betweenness',
                              ifelse((impact < quantile(impact,.75)) & (log(btw+1) > quantile(log(btw+1),.75)),'Low impact, high betweenness',
                                     ifelse((impact > quantile(impact,.75)) & (log(btw+1) < quantile(log(btw+1),.75)),'High impact, low betweenness','High impact, high betweenness')))))

  rownames(df.res) = data$Knockout[index]
  q1 <- row.names(subset(df.res, (impact < quantile(impact,.75)) & (log.betweenness < quantile(log.betweenness,.75))))
  q2 <- row.names(subset(df.res, (impact < quantile(impact,.75)) & (log.betweenness > quantile(log.betweenness,.75))))
  q3 <- row.names(subset(df.res, (impact > quantile(impact,.75)) & (log.betweenness < quantile(log.betweenness,.75))))
  q4 <- row.names(subset(df.res, (impact > quantile(impact,.75)) & (log.betweenness > quantile(log.betweenness,.75))))
  idx <- unique(c(sample(q1,6),sample(q2,6),sample(q3,6),sample(q4,6)))
  df.idx <- df.res[idx,]

  gp1 <-  ggplot2::ggplot(data=df.res,ggplot2::aes(x=impact,y=log.betweenness)) +
    geom_point(ggplot2::aes(col=pos.label), alpha=.3, size=3) +
    scale_color_manual(values=c("plum4","palegreen4","indianred","cornflowerblue")) +
    theme_linedraw() +
    theme_light() +
    theme(legend.position="bottom") +
    guides(colour = guide_legend(nrow = 2)) +
    theme(legend.title=element_blank()) +
    geom_text_repel(data=df.idx, ggplot2::aes(label=df.idx$Knockout), size=3.5) +
    geom_vline(xintercept=quantile(df.res$impact,.75),linetype="dashed") +
    geom_hline(yintercept=quantile(df.res$log.betweenness,.75),linetype="dashed") +
    xlab("Impact") +
    ylab("Log(betweenness+1)")


  rownames(df.res) <- c()
  df.res2 <- df.res[,-c(4,5)]
  names(df.res2) <- c('Knockout','Impact','Betweenness','Position')

  gene.cluster <- df.symb[,c(1,18)]
  names(gene.cluster) <- c('Knockout','Cluster')
  df.res3 <-merge(df.res2, gene.cluster, by="Knockout",all.x=TRUE)

  df.tab <- data.frame(table(df.res3$Cluster,df.res3$Position))
  df.tab2 <- df.tab %>% group_by(Var1) %>% top_n(1, Freq)
  names(df.tab2) <- c('Cluster','Position','nGenes')


  #### -------------------> Output
  res <- list()
  class(res) = "GeneNetwork"
  res$plot.pnet <- gNEt(dat=Net) # plot gene network
  res$plot.impact_betweenees <- gp1 # plot impact betweenees
  res$stats.impact_betweeness <- df.res3 # impact betweenees data
  res$stats.impact_betweeness_by_cluster <- df.tab2 # plot position by cluster

  return(res)
}




