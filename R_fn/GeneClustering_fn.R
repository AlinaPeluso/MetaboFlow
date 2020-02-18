
Packages <- c("dplyr","tidyr","reshape2","reshape","ggplot2","ggfortify","factoextra","pheatmap","gplots","mixOmics","qgraph","GOstats","GO.db","org.Sc.sgd.db","data.table","gridExtra","sna","GGally","intergraph","igraph","network","Matrix","ggrepel","knitr","tidyr","psych")
invisible(lapply(Packages, library, character.only = TRUE))

GeneClustering = function(data=NULL,
                          data_Symb=NULL)
{

  data_GOslim <- read.table('./go_slim_mapping.txt', sep='\t', header=T)
  data_ORF2KEGG <- read.table('./ORF2KEGG.txt', sep='\t', header=T)
  

  #### -------------------> Define clusters
  res.dist <- dist(data_Symb[,-1], method = "manhattan")
  res.hc <- hclust(d = res.dist, method = "single")
  data_Symb$cluster <- cutree(res.hc, h = 0) # distance 0
  
  #### -------------------> Subset cluster with more than 10 genes
  df <- as.data.frame(table(data_Symb$cluster)); names(df) <- c('cluster', 'nGenes')
  df_sub <- df[df$nGenes>10,]
  rownames(df_sub) <- c()
  
  #### -------------------> Plot clustering
  df.tmp_long = list()
  for(i in 1:dim(df_sub)[1]){
    tmp_wide <- data[data_Symb$cluster==df_sub$cluster[i],]
    tmp_long <- gather(tmp_wide, Ion, logConcentration_corr_norm,  Ca:Zn, factor_key=TRUE)
    tmp_long$Cluster_number <- rep(df_sub$cluster[i],nrow(tmp_long))
    label <- paste("Cluster",df_sub$cluster[i],paste("(",df_sub$nGenes[i], " genes)", sep=''), sep=' ')
    tmp_long$Cluster <- rep(label,nrow(tmp_long))
    df.tmp_long[[i]] <- tmp_long
  } 
  df.tmp_long_all <- data.frame(do.call(rbind,df.tmp_long))
  
  p1 <- ggplot(data = df.tmp_long_all, aes(x = Ion, y = logConcentration_corr_norm)) +
    facet_wrap(~Cluster) + 
    geom_line(aes(group = Knockout)) +
    stat_summary(fun.data = "mean_se", color = "red")+
    labs(x = "", y = "") +
    theme(legend.position="none") +
    theme(axis.text.x=element_text(angle=90, hjust=1),axis.text=element_text(size=10))
  
  
  
  #### -------------------> KEGG AND GO SLIM ANNOTATION
  p.data_list = list()
  for(i in 1:dim(df_sub)[1]){
    inputGeneSet <- data_Symb$Knockout[data_Symb$cluster==df_sub$cluster[i]]
    
    N <- as.numeric(length(inputGeneSet))
    p.data <<- data_GOslim %>%
      dplyr::mutate(Ontology = setNames(c("Biological process","Cellular component","Molecular function"),c("P","C","F"))[as.character(Ontology)]) %>%
      dplyr::filter(ORFs %in% as.character(inputGeneSet)) %>%
      dplyr::group_by(GOslim,Ontology) %>%
      dplyr::filter(GOslim != "other") %>%
      dplyr::rename(Term = GOslim) %>%  
      dplyr::summarise(Count = n()) %>%
      dplyr::mutate(Percent = Count/N*100) %>%
      dplyr::bind_rows(data_ORF2KEGG %>%
                         dplyr::filter(ORF %in% as.character(inputGeneSet)) %>%
                         dplyr::group_by(KEGGID, Pathway) %>%
                         dplyr::summarise(Count = dplyr::n()) %>%
                         dplyr::mutate(Ontology = "KEGG") %>%
                         dplyr::rename(Term = Pathway) %>%
                         dplyr::ungroup() %>%
                         dplyr::filter(KEGGID != "01100") %>%
                         dplyr::select(-KEGGID)%>%
                         dplyr::mutate(Percent = Count/N*100)) %>%
      dplyr::filter(!Term %in% c("molecular_function", "biological_process", "cellular_component"))
    p.data_list[[i]] = p.data
  }
  
  
  df.data_list = label_list <- list()
  for(i in 1:dim(df_sub)[1]){
    p.data_list_sub <- as.data.frame(p.data_list[[i]])
    p.data_list_sub$Cluster_number <- rep(df_sub$cluster[i],nrow(p.data_list_sub))
    label <- paste("Cluster",df_sub$cluster[i],paste("(",df_sub$nGenes[i], " genes)", sep=''), sep=' ')
    p.data_list_sub$Cluster <- rep(label,nrow(p.data_list_sub))
    sub_data <- p.data_list_sub[(p.data_list_sub$Percent>5) & (p.data_list_sub$Ontology %in% c("Biological process","Cellular component","Molecular function")),]
    df.data_list[[i]] <- sub_data
    label_list[[i]] <- label
  }  
  
  Kegg_Goslim_annotation.list_by_clust <- lapply(df.data_list,"[", c('Term','Ontology','Count','Percent'))
  names(Kegg_Goslim_annotation.list_by_clust) <- label_list
  
  #### -------------------> GO TERMS ENRICHMENT
  universeGenes=as.character(data_Symb$Knockout) # universo dei dati sperimentali
  
  p.data_list2 = label_list2 = list()
  for(i in 1:dim(df_sub)[1]){
    label <- paste("Cluster",df_sub$cluster[i],paste("(",df_sub$nGenes[i], " genes)", sep=''), sep=' ')
    inputGeneSet <- data_Symb$Knockout[data_Symb$cluster==df_sub$cluster[i]] 
    
    ont=c("BP","MF","CC")
    results<<-c()
    for(k in 1:3){
      params = new("GOHyperGParams",
                   geneIds=as.character(inputGeneSet),
                   universeGeneIds=universeGenes,
                   annotation="org.Sc.sgd.db",
                   categoryName="GO",
                   ontology=ont[k],
                   pvalueCutoff=0.05,
                   conditional=T,
                   testDirection="over")
      hgOver <<- hyperGTest(params)
    }  
    results<<-rbind(results,cbind(setNames(dplyr::data_frame(ID=names(pvalues(hgOver)), Term = Term(ID),
                                                             pvalues = pvalues(hgOver), oddsRatios = oddsRatios(hgOver),expectedCounts=expectedCounts(hgOver), 
                                                             geneCounts=geneCounts(hgOver),universeCounts=universeCounts(hgOver)),
                                           c("GO_ID","Description","Pvalue","OddsRatio","ExpCount","Count","CountUniverse")),"Ontology"=ont[k])
                    %>%dplyr::filter(Pvalue <= 0.05&Count > 1))
    
    p.data_list2[[i]] <- results
    label_list2[[i]] <- label
  }
  Goterms_enrichment.list_by_clust <- lapply(p.data_list2,"[",-c(4,5,8))
  names(Goterms_enrichment.list_by_clust) <- label_list2
  
  
  #### -------------------> Output
  res <- list()
  class(res) = "GeneClustering"
  
  res$stats.clusters <- df_sub # selected clusters (nGenes>10)
  res$plot.profiles <- p1 # plot cluster profiles
  res$stat.Kegg_Goslim_annotation <- Kegg_Goslim_annotation.list_by_clust # KEGG AND GO SLIM ANNOTATION
  res$stat.Goterms_enrichment <- Goterms_enrichment.list_by_clust # GO TERMS ENRICHMENT
  
  return(res)
}




