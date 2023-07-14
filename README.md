# cNFOrganoidAnalysis
This repository is designed to analyze data from a panel of [cNF Organoids](https://www.synapse.org/#!Synapse:syn11374354/wiki/488832). We are currently evaluating the data to assess the fidelity off the organoids to the primary patient data. All of the data is stored on Synapse, which is the back end to the NF Data Portal. This repository serves two purposes: first, to store our work as we analyze and assess which results to present, and second, to enable reproducibility of our figures. 

Here we focus on the files that enable the figure, specifically Figures 4 and 5 of our manuscript, to be recreated. 

## Helper files
Our helper files do not create figures directly but carry out some basic functions to standardize the data retreival and ploting functions.

|File| Description|
|---|---|
|[./loadExpFromCounts.R](loadExpFromCounts.R)| This is an R script that loads data from Synaspse|
|[./orgPlottingFunctions.R](orgPlottingFunctions.R)| THis contains plotting functions|


## Figures
Figures 4 and 5 of the manuscript were generated from the sequencing data, and therefore have scripts to generate them.

|File| Description|
|---|---|
|[./fig4ab_plotGeneExpression.R](fig4ab_plotGeneExpression.R)|Plots panels A and B of figure 4|
|[./fig4cd_tumorDeconv.R](fig4cd_tumorDeconv.R)| Tjis contains plotting functions|
|[./fig5_meta_analysis.R](fig5_meta_analysis.R)| This is the analysis comparing correlations across data modalities.|

## Samples
We combine samples from the CTF Biobank as well as our own, depicted below.
[!samples](./pcaOfAllSamples.png)

