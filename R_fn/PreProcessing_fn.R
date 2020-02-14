
Packages <- c("dplyr","tidyr","reshape2","reshape","ggplot2","ggfortify","factoextra","pheatmap","gplots","data.table","gridExtra","network","Matrix","ggrepel","knitr","tidyr","psych")
invisible(lapply(Packages, library, character.only = TRUE))

PreProcessing = function(data=NULL,stdev=NULL)
{

  #### -------------------> Import data 
  list.s <- list()
  data2 <- data[,-c(1,2)]
  for (i in 1:14){
    list.s[[i]] <- c(names(data2)[i],round(summary(data2[,i]),3),round(var(data2[,i]),3))
  }
  df.s <- data.frame(do.call(rbind,list.s))
  names(df.s) <-  c('Ion','Min','1st Quartile','Median','Mean', '3rd Quartile', 'Max','Variance' )
  
  #### -------------------> Outlier detection
  data$id <- row.names(data)
  data_long <- gather(data,Ion,Concentration, Ca:Zn, factor_key=TRUE) 
  
  for (i in 1:length(levels(data_long$Ion))){
    data_long_sub <- data_long[data_long$Ion==levels(data_long$Ion)[i],]
    lowerq = quantile(data_long_sub$Concentration,na.rm=T)[2] 
    upperq = quantile(data_long_sub$Concentration,na.rm=T)[4]
    iqr = upperq - lowerq 
    extreme.t.upper = (iqr * 3) + upperq
    extreme.t.lower = lowerq - (iqr * 3)
    data_long[data_long$Ion==levels(data_long$Ion)[i],'Outlier'] <- ifelse((data_long_sub$Concentration > extreme.t.upper | data_long_sub$Concentration < extreme.t.lower),1,0)
  }
  df_outlier <- data.frame(cbind(levels(data_long$Ion),table(data_long$Ion,data_long$Outlier),round(table(data_long$Ion,data_long$Outlier)[,2]/dim(data_long)[1]*100,2)))
  rownames(df_outlier) <- c()
  colnames(df_outlier) <- c('Ion','no outlier','outlier','outlier(%)')
  data_long_clean <- data_long[data_long$Outlier<1,] 
  
  #### -------------------> Median batch correction
  data_long_clean$logConcentration <- log(data_long_clean$Concentration)
  med <- matrix(0,length(table(data$Batch_ID)),length(levels(data_long$Ion)))
  datasets <- list()
  j <- 1
  while (j <= length(levels(data_long_clean$Ion))){
    data_long_clean_sub <- data_long_clean[data_long_clean$Ion==levels(data_long_clean$Ion)[j],]
    for (i in 1:max(data_long_clean_sub$Batch_ID)){
      med[i,j] <- median(data_long_clean_sub$logConcentration[data_long_clean_sub$Batch_ID==i])
      data_long_clean_sub$logConcentration_corr[data_long_clean_sub$Batch_ID==i] <- data_long_clean_sub$logConcentration[data_long_clean_sub$Batch_ID==i] - med[i,j]
    }
    datasets[[j]] <- data_long_clean_sub
    data_long_clean_scaled <- do.call(rbind, datasets)
    j <- j+1
  }
  
  p1 <- ggplot(data = data_long_clean_scaled, aes(x = factor(Batch_ID), y = logConcentration_corr, col=factor(Batch_ID))) + 
    geom_point(shape=1) +
    facet_wrap(~Ion) + 
    xlab("Batch.ID") + 
    ylab("log(Concentration) (ppm)") +
    theme(legend.position="none") +
    theme(axis.text.x=element_blank())
  
  
  list.mbc <- list()
  for (i in 1:max(length(levels(data_long_clean_scaled$Ion)))){
    list.mbc[[i]] <- c(levels(data_long_clean_scaled$Ion)[i],round(summary(datasets[[i]]$logConcentration_corr),3),round(var(datasets[[i]]$logConcentration_corr),3))
  }
  df.mbc <- data.frame(do.call(rbind,list.mbc))
  names(df.mbc) <- c('Ion','Min','1st Quartile','Median','Mean', '3rd Quartile', 'Max','Variance' )
  
  #### -------------------> Standardisation
  if (is.null(stdev)) {
    sds <- rep(NA,length(levels(data_long_clean_scaled$Ion)))
    for (i in 1:length(levels(data_long_clean_scaled$Ion))){
      sds[i] <- as.numeric(sd(data_long_clean_scaled$logConcentration_corr[data_long_clean_scaled$Ion==levels(data_long_clean_scaled$Ion)[i]])) # Ions' sd of logConcentration_corr 
    }
  }else{
    sds=as.numeric(as.vector(stdev[,1]))

  }
  
  datasets2 <- list()
  j <- 1
  while (j <= length(levels(data_long_clean_scaled$Knockout))){     
    data_long_clean_scaled_sub <- data_long_clean_scaled[data_long_clean_scaled$Knockout==levels(data_long_clean_scaled$Knockout)[j],]
    for (i in 1:length(levels(data_long_clean_scaled$Ion))){
      data_long_clean_scaled_sub$logConcentration_corr_norm[data_long_clean_scaled_sub$Ion==levels(data_long_clean_scaled_sub$Ion)[i]] <- 
        data_long_clean_scaled_sub$logConcentration_corr[data_long_clean_scaled_sub$Ion==levels(data_long_clean_scaled_sub$Ion)[i]] / sds[i]
    }
    datasets2[[j]] <- data_long_clean_scaled_sub
    data_long_clean_scaled_norm <- do.call(rbind, datasets2)
    j <- j+1
  }
  
  list.mbc2 <- list()
  for (i in 1:max(length(levels(data_long_clean_scaled_norm$Ion)))){
    list.mbc2[[i]] <- c(levels(data_long_clean_scaled_norm$Ion)[i],round(summary(datasets2[[i]]$logConcentration_corr_norm),3),round(var(datasets2[[i]]$logConcentration_corr_norm),3))
  }
  
  df.mbc2 <- data.frame(do.call(rbind,list.mbc2))
  names(df.mbc2) <- c('Ion','Min','1st Quartile','Median','Mean', '3rd Quartile', 'Max','Variance' )
  
  
  #### -------------------> Symbolization
  data_long_clean_scaled_norm$Symb <- ifelse((data_long_clean_scaled_norm$logConcentration_corr_norm > -3) & (data_long_clean_scaled_norm$logConcentration_corr_norm< 3), 0, ifelse(data_long_clean_scaled_norm$logConcentration_corr_norm>=3,1,-1)) 
  
  #### -------------------> Aggregation of the batch replicas
  data_long_clean_scaled_norm_unique <- data.frame(aggregate(. ~ Knockout*Ion, data_long_clean_scaled_norm[,c('Knockout','Ion','logConcentration_corr_norm','Symb')], median))
  data_long_clean_scaled_norm_unique$Symb <- ifelse((data_long_clean_scaled_norm_unique$Symb<0.5) & (data_long_clean_scaled_norm_unique$Symb>-0.5), 0, ifelse(data_long_clean_scaled_norm_unique$Symb>=0.5,1,-1)) 
  data_wide_clean_scaled_norm_unique <- reshape2::dcast(data_long_clean_scaled_norm_unique, Knockout~ Ion, value.var="logConcentration_corr_norm")
  data_wide_clean_scaled_norm_unique_Symb <- reshape2::dcast(data_long_clean_scaled_norm_unique, Knockout~ Ion, value.var="Symb")
  
  p2 <- ggplot(data = data_long_clean_scaled_norm_unique, aes(x = logConcentration_corr_norm)) + 
    geom_histogram(binwidth=.1) +
    facet_wrap(~Ion) + 
    xlab("log(Concentration) (z-score)") + 
    ylab("Frequency") + 
    geom_vline(xintercept=c(-3,3),col='red') 
  
  
  
  #### -------------------> Output
  res <- list()
  class(res) = "PreProcessing"
  
  res$stats.raw_data <- df.s # raw data
  res$stats.outliers <- df_outlier # outliers
  res$stats.median_batch_corrected_data <- df.mbc # median batch corrected data
  res$stats.standardised_data <- df.mbc2 # standardised data
  
  res$dataR.long <- data_long_clean_scaled_norm
  res$data.long <- data_long_clean_scaled_norm_unique
  res$data.wide <- data_wide_clean_scaled_norm_unique
  res$data.wide_Symb <- data_wide_clean_scaled_norm_unique_Symb
  
  res$plot.logConcentration_by_batch <- p1
  res$plot.logConcentration_z_scores <- p2
  return(res)
}

