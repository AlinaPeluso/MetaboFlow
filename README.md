# MetaboFlow
This pipeline describes the processing and the analysis of Ionomics data. 
This [paper](https://arxiv.org/abs/1910.14191) describes a possible application.

The workflow consist of four sections, and respectively

* Pre-processing
* Exploratory analysis
* Clustering which also includes the GO Slim annotation and the GO terms enrichment
* Network analysis

The pre-processing section is required as first as it produces the cleaned dataset used in the other sections. There is no specific order on how to run the other three sections. 

## Pre-processing
This first section aims to free the data from unreliable samples which will probably lead to wrong outputs. In such way, effective data pre-processing methods are applied to avoid the effects of noisy and unreliable data.

This section requires as input the raw data frame, e.g. ion's concentrations. It is also possible to define a set of ion's standard deviation, as these are possibly computed accounting for some control genes. Note that the latter is an optional input i.e if not provided the standard deviations from the data would be computed to perform the data standardisation (see Section [Standardisation](https://github.com/AlinaPeluso/MetaboFlow#standardisation).


```
# Import data
Idata <- read.csv('./Data/Dataset_Metaboflow_Ionomic_Workflow.csv',header=T)
pre_defined_sd <- read.table('./Data/pre_defined_sd.txt', header=T)

# Execute pre-processing fn
source("./R_fn/PreProcessing_fn.R")
fn1 <- PreProcessing(data=Idata,stdev=pre_defined_sd)
```

#### Inspect the imported raw data

We illustrate the Ionomics workflow with ICP-MS data of yeast intracellular ion concentrations measured for 1454 single-gene haploid knockouts [see @danku2009high]. Ions measured include Ca44, Cd111, Co59, Cu65, Fe56, K39, Mg25, Mn55, Mo95, Na23, Ni60, P31, S34, Zn66. Values of concentration are in ppm and have been adjusted using optical density measurements. Intracellular concentrations are measured for two, four or eight replicas of each mutant depending on the knock-out. Knock-out YDL227C mutant is measured multiple times in every batch as control strain.

The concentration values for the raw data ion can be summarised as follow.

```
fn1$stats.raw_data
```

|    | Ion | Min       | 1st Quartile | Median    | Mean      | 3rd Quartile | Max        | Variance    |
|----|-----|-----------|--------------|-----------|-----------|--------------|------------|-------------|
| 1  | Ca  | 0\.449    | 31\.73       | 40\.44    | 45\.071   | 51\.015      | 902\.568   | 829\.525    |
| 2  | Cd  | 0\.174    | 0\.866       | 0\.988    | 1\.002    | 1\.121       | 2\.512     | 0\.051      |
| 3  | Co  | 0\.007    | 0\.142       | 0\.16     | 0\.162    | 0\.184       | 0\.702     | 0\.001      |
| 4  | Cu  | 0\.587    | 1\.344       | 1\.586    | 1\.717    | 1\.831       | 327\.79    | 16\.91      |
| 5  | Fe  | 0\.002    | 5\.527       | 7\.295    | 9\.469    | 9\.332       | 6624\.526  | 5154\.611   |
| 6  | K   | 284\.273  | 2060\.619    | 2495\.265 | 2492\.765 | 2879\.551    | 17777\.452 | 534784\.375 |
| 7  | Mg  | 115\.63   | 546\.325     | 679\.275  | 642\.578  | 753\.678     | 3838\.479  | 31947\.598  |
| 8  | Mn  | 0\.02     | 0\.982       | 1\.206    | 1\.197    | 1\.38        | 7\.339     | 0\.106      |
| 9  | Mo  | 0\.158    | 0\.656       | 0\.934    | 1\.109    | 1\.327       | 60\.879    | 1\.855      |
| 10 | Na  | 0\.184    | 128\.831     | 185\.25   | 196\.545  | 247\.747     | 892\.968   | 9027\.944   |
| 11 | Ni  | 0\.074    | 0\.982       | 1\.258    | 1\.693    | 1\.543       | 2323\.058  | 565\.618    |
| 12 | P   | 1194\.953 | 3833\.9      | 4514\.476 | 4289\.109 | 4952\.98     | 21695\.748 | 1151197\.05 |
| 13 | S   | 20\.592   | 434\.61      | 512\.845  | 529\.493  | 605\.436     | 5484\.638  | 37212\.137  |
| 14 | Zn  | 7\.659    | 14\.785      | 16\.549   | 17\.114   | 18\.334      | 2221\.586  | 511\.804    |


There is a very high variability of the knockouts across ions and within the batches.
There are no missing values in the data.


#### Outlier detection

We define a lower outer fence: `$Q1 - 3*IQ$` and a upper outer fence: `$Q3 + 3*IQ$` where `$Q1$` and `$Q3$` are the first and the third quantile of the distribution, respectively. A point beyond the outer fence is considered an extreme outlier.

The outliers are split across ions as follows. 

```
fn1$stats.outliers
```
| Ion | no | outlier | outlier | outlier\(%\) |
|-----|----|---------|---------|--------------|
| 1   | Ca | 9694    | 305     | 0\.22        |
| 2   | Cd | 9950    | 49      | 0\.04        |
| 3   | Co | 9966    | 33      | 0\.02        |
| 4   | Cu | 9870    | 129     | 0\.09        |
| 5   | Fe | 9833    | 166     | 0\.12        |
| 6   | K  | 9980    | 19      | 0\.01        |
| 7   | Mg | 9991    | 8       | 0\.01        |
| 8   | Mn | 9984    | 15      | 0\.01        |
| 9   | Mo | 9909    | 90      | 0\.06        |
| 10  | Na | 9965    | 34      | 0\.02        |
| 11  | Ni | 9849    | 150     | 0\.11        |
| 12  | P  | 9991    | 8       | 0\.01        |
| 13  | S  | 9923    | 76      | 0\.05        |
| 14  | Zn | 9888    | 111     | 0\.08        |


#### Median batch correction

First we take the logarithmm of the concentration value. Then, the data are scaled to the median taken for each ion within each batch.

```
fn1$stats.median_batch_corrected_data
```

After outlier removal and the median batch correction of the logged concentrations (logConcentration_corr), the data looks as

```
fn1$plot.logConcentration_by_batch 
```
![caption](link)

#### Standardisation
