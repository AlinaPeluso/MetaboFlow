
install.packages("devtools") 
library(devtools)
install_github("AlinaPeluso/MetaboFlow", subdir="/IonFlow")
library(IonFlow)

### ***Pre-processing***
pre_proc <- PreProcessing(data=IonData,stdev=pre_defined_sd)
# stats
pre_proc$stats.raw_data 
pre_proc$stats.outliers 
pre_proc$stats.median_batch_corrected_data
pre_proc$stats.standardised_data 
# plots
pre_proc$plot.logConcentration_by_batch
pre_proc$plot.logConcentration_z_scores
# data
head(pre_proc$dataR.long)
head(pre_proc$data.long)
head(pre_proc$data.wide)
head(pre_proc$data.wide_Symb)

### ***Exploratory analysis***
exp_anal <- ExploratoryAnalysis(data=pre_proc$data.wide)
# plots
exp_anal$plot.Pearson_correlation
exp_anal$plot.PCA_Individual
exp_anal$plot.heatmap
exp_anal$plot.pairwise_correlation_map
exp_anal$plot.regularized_partial_correlation_network
# data
head(exp_anal$data.PCA_loadings)

### ***Gene Clustering***
gene_clust <- GeneClustering(data=pre_proc$data.wide, data_Symb=pre_proc$data.wide_Symb)
# stats
gene_clust$stats.clusters
gene_clust$stats.Kegg_Goslim_annotation
gene_clust$stats.Goterms_enrichment
# plots
gene_clust$plot.profiles

### ***Gene Network***
gene_net <- GeneNetwork(data=pre_proc$data.wide, data_Symb=pre_proc$data.wide_Symb)
# stats
gene_net$stats.impact_betweeness
gene_net$stats.impact_betweeness_by_cluster
# plots
gene_net$plot.pnet 
gene_net$plot.impact_betweenees
