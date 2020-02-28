---
title: 'IonFlow: Pipeline for processing and analysis of Ionomics data'
tags:
  - R
  - Ionomics
  - Metabolomics
  - MetaboFlow
authors:
  - name: Alina Peluso
    orcid: 0000-0003-2895-0406
    affiliation: "1, 2" 
  - name: Jacopo Iacovacci
    affiliation: "1, 2"
  - name: Robert Glen
    affiliation: "1, 3"
affiliations:
 - name: Department of Metabolism Digestion and Reproduction, Faculty of Medicine, Imperial College London, London, SW7 2AZ, United Kingdom
   index: 1
 - name: The Molecular Biology of Metabolism Laboratory,
The Francis Crick Institute, London, NW1 1AT, United Kingdom
   index: 2
 - name: Centre for Molecular Informatics, Department of Chemistry,
University of Cambridge, Lensfield Road, Cambridge, CB21EW, United Kingdom
   index: 3
date: 27 February 2020
bibliography: paper.bib
---

# Summary

This package illustrate tools for processing Ionomics data to aid reproducible data sharing and big data initiatives. 

The workflow is available within the `IonFlow` R package, and it consists of four sections:

* Pre-processing,
* Exploratory analysis,
* Clustering which also includes the GO Slim annotation and GO terms enrichment,
* Network analysis.

To show the package capabilities we illustrate the Ionomics workflow employing ICP-MS data of yeast intracellular ion concentrations measured for a set of single-gene haploid knockouts [see @danku2009high]. Ions measured include Ca44, Cd111, Co59, Cu65, Fe56, K39, Mg25, Mn55, Mo95, Na23, Ni60, P31, S34, Zn66. Values of concentration are in ppm and have been adjusted using optical density measurements. Intracellular concentrations are measured for two, four or eight replicas of each mutant depending on the knock-out. Knock-out YDL227C mutant is measured multiple times in every batch as control strain.

#### Pre-processing 
The pre-processing section is required first as it produces in output the cleaned dataset to be used in the other sections. The aim of this first section is to free the data from unreliable samples which will probably lead to wrong outputs. In such way, effective data pre-processing methods are applied to avoid the effects of noisy and unreliable data. This section includes outliers detection and removal, median batch correction, standardisation, symbolisation of the ionomics profiles, and aggregation of the genes knockout replicas. Three dataset are obtained as output. The first in the long format (genes knockout as rows and ions as columns), and two in wide format and respectively one with the standardised ion's concentraction, and the other with the symbolised profiles of the knockouts.

#### Exploratory analysis

The exploratory analysis section provides a way to summarize the main characteristics of the data with visual methods. This full section can be run employong the datasets obtained as output by the pre-processing section. This section includes profiling of (partial) correlations of interest, dimensionality reduction methods such as principal component analysis, and graphical representations of the data such as heatmap of the genes knockout with a clustering dendrogram.

#### Clustering

The clustering algorithm is based on the computation of the Manhattan distance between the knockouts' symbolised profile to cluster genes having same symbolic profile, that is  the unique partition is found by cutting the hierarchical tree at zero-distances. A profile plot is generated to investigate the components of each cluster. Next, we employ the Go Slim annotation to highlight the biological process, the cellular component, and the molecular function of the genes within each cluster. In particular, we retain the annotations that map at least 5% of the genes in the cluster. 
We also perform the Go terms enrichment of the genes in the clusters by employing the GO terms annotation in the [SGD online database](https://www.yeastgenome.org/) [see @engel2014reference].


#### Network analysis

The network analysis is based on the computation of the empirical correlation matrix between genes. Moreover, as we are interested only in positive correlations among genes we filter the correlation matrix based on an arbitrary correlation's thresold of 0.6 such that the clustering can be interpreted in terms of posivite correlation between gene profiles.
A network plot is generated to visualise the network structure. Next, we 
focus on the identification of the most central genes within each network subspace. To do so we can consider two metrics i.e. the impact and the betweeness score. 
From the empirical correlation matrix between genes we can compute the betweenness measure as the fraction of shortest paths that pass through each gene (node). Next, a measure of impact can be computed as the $L_2$ norm (Euclidean distance) of each gene. These two centrality measures can be then used togheter to cluster the genes. In addition to this, we can also associate each cluster to low or high values of impact and betwenees based on the highest number of genes in that cathegory. 



# Acknowledgment

The authors acknowledge funding from the Wellcome Trust funded project MetaboFlow, grant reference number 202952/D/16/Z.