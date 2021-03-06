% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/GeneNetwork_fn.R
\name{GeneNetwork}
\alias{GeneNetwork}
\title{Network}
\usage{
GeneNetwork(data = NULL, data_Symb = NULL)
}
\arguments{
\item{\code{data}}{a processed data frame of ion's concentrations (columns) for each gene's knockout (rows) (wide format) e.g. \code{data.wide} from the \code{PreProcessing_fn}}

\item{\code{data_Symb}}{a processed data frame of symbolised (values -1,0, or 1) ion's concentration profiles (columns) for each gene's knockout (rows) (wide format) e.g. \code{data.wide_Symb} from the \code{PreProcessing_fn}}
}
\value{
\code{stats.impact_betweeness} statistics of the impact betweenees data

\code{stats.impact_betweeness_by_cluster} statistics of the impact betweenees data by cluster

\code{plot.pnet plot} of the gene's network

\code{plot.impact_betweenees} plot of the impact betweenees
}
\description{
Network tools for Ionomics data
}
\examples{
\code{
library(IonFlow)

### Run PreProcessing function
pre_proc <- PreProcessing(data=IonData,stdev=pre_defined_sd)

### Run GeneNetwork function
gene_net <- GeneNetwork(data=pre_proc$data.wide, data_Symb=pre_proc$data.wide_Symb)

# stats
gene_net$stats.impact_betweeness
gene_net$stats.impact_betweeness_by_cluster

# plots
gene_net$plot.pnet
gene_net$plot.impact_betweenees
}
}
