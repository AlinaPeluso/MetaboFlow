

#====================================================================================================
#### f1: PreProcessing

# *** Inputs ***
Idata <- read.csv('./Data/Dataset_Metaboflow_Ionomic_Workflow.csv',header=T)
pre_defined_sd <- read.table('./Data/pre_defined_sd.txt', header=T)

# *** Run function ***
f1 <- PreProcessing(data=Idata,stdev=pre_defined_sd)


# *** Outputs ***
f1$stats.raw_data
f1$stats.outliers
f1$stats.median_batch_corrected_data
f1$stats.standardised_data

head(f1$dataR.long)
head(f1$data.long)
head(f1$data.wide)
head(f1$data.wide_Symb)

f1$plot.logConcentration_by_batch 
f1$plot.logConcentration_z_scores



# *** Store outputs ***
setwd("./Output/f1/")

cat(capture.output(print(f1$stats.raw_data), file="stats.raw_data.txt"))
cat(capture.output(print(f1$stats.outliers), file="stats.outliers.txt"))
cat(capture.output(print(f1$stats.median_batch_corrected_data), file="stats.median_batch_corrected_data.txt"))
cat(capture.output(print(f1$stats.standardised_data), file="stats.standardised_data.txt"))

write.csv(f1$data.long, file = "data.long.csv", row.names = T)
write.csv(f1$data.wide, file = "data.wide.csv", row.names = T)
write.csv(f1$data.wide_Symb, file = "data.wide_Symb.csv", row.names = T)

pdf(file = "plot.logConcentration_by_batch.pdf",) 
f1$plot.logConcentration_by_batch 
dev.off()

pdf(file = "plot.logConcentration_z_scores.pdf",) 
f1$plot.logConcentration_z_scores
dev.off()


#====================================================================================================
#### f2: ExploratoryAnalysis


# *** Inputs ***
#data.wide <- read.table('./Idata_wide_clean_scaled_norm_unique.txt', sep='\t', header=T)
#f2 <- ExploratoryAnalysis(data=data.wide) 


# *** Run function ***
f2 <- ExploratoryAnalysis(data=f1$data.wide) 

# *** Outputs ***
f2$plot.Pearson_correlation 
f2$plot.PCA_Individual

f2$data.PCA_loadings

f2$plot.heatmap 
f2$plot.pairwise_correlation_map 
f2$plot.regularized_partial_correlation_network 


# *** Store outputs ***
setwd("./Output/f2/")

pdf(file = "plot.Pearson_correlation.pdf",) 
f2$plot.Pearson_correlation 
dev.off()

pdf(file = "plot.PCA_Individual.pdf",) 
f2$plot.PCA_Individual
dev.off()

write.csv(f2$data.PCA_loadings , file = "data.PCA_loadings.csv", row.names = T)

pdf(file = "plot.heatmap.pdf",) 
f2$plot.heatmap 
dev.off()

pdf(file = "plot.pairwise_correlation_map.pdf",) 
f2$plot.pairwise_correlation_map
dev.off()

pdf(file = "plot.regularized_partial_correlation_network.pdf",) 
f2$plot.regularized_partial_correlation_network 
dev.off()



#====================================================================================================
#### f3: GeneClustering



# *** Inputs ***
go_slim_mapping <- read.table('./go_slim_mapping.txt', sep='\t', header=T)
ORF2KEGG <- read.table('./ORF2KEGG.txt', sep='\t', header=T)
#data.wide <- read.table('./Idata_wide_clean_scaled_norm_unique.txt', sep='\t', header=T)
#data.wide_Symb <- read.table('./Idata_wide_clean_scaled_norm_unique_Symb.txt', sep='\t', header=T)
#f3 <- GeneClustering(data=data.wide, data_Symb=data.wide_Symb) 

# *** Run function ***
f3 <- GeneClustering(data=f1$data.wide, data_Symb=f1$data.wide_Symb) 


# *** Outputs ***
f3$stats.clusters 

f3$plot.profiles 

f3$stat.Kegg_Goslim_annotation 
f3$stat.Goterms_enrichment 

# *** Store outputs ***
setwd("./Output/f3/")

cat(capture.output(print(f3$stats.clusters), file="stats.clusters.txt"))

pdf(file = "plot.profiles.pdf",) 
f3$plot.profiles 
dev.off()

cat(capture.output(print(f3$stat.Kegg_Goslim_annotation), file="stat.Kegg_Goslim_annotation.txt"))
cat(capture.output(print(f3$stat.Goterms_enrichment), file="stat.Goterms_enrichment.txt"))

#====================================================================================================
#### f4: GeneNetwork


# *** Inputs ***
#data.wide <- read.table('./Idata_wide_clean_scaled_norm_unique.txt', sep='\t', header=T)
#data.wide_Symb <- read.table('./Idata_wide_clean_scaled_norm_unique_Symb.txt', sep='\t', header=T)
#f4 <- GeneNetwork(data=data.wide, data_Symb=data.wide_Symb) 

# *** Run function ***
f4 <- GeneNetwork(data=f1$data.wide, data_Symb=f1$data.wide_Symb) 


# *** Outputs ***
f4$plot.pnet
f4$plot.impact_betweenees 

f4$stats.impact_betweeness 
f4$stats.impact_betweeness_by_cluster


# *** Store outputs ***
setwd("./Output/f4/")

pdf(file = "plot.pnet.pdf",) 
f4$plot.pnet
dev.off()

pdf(file = "plot.impact_betweenees .pdf",) 
f4$plot.impact_betweenees 
dev.off()

cat(capture.output(print(f4$stats.impact_betweeness), file="stats.impact_betweeness.txt"))
cat(capture.output(print(f4$stats.impact_betweeness_by_cluster), file="stats.impact_betweeness_by_cluster.txt"))


