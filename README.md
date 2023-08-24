---
title: "MRS"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

# MRS

Welcome to *MRS*'s page!


Installing *MRS*
-------------
To install the development version from github:

```r
library(devtools)
install_github("chanw0/MRS")
```
Once installed, the package can be loaded using:

```r
library("MRS")
```

## Examples


### Example 1: Evaluation of MRS in terms of comparison between Healthy and Nonhealthy 

1.1 using ANCOMBC method and Shannon index 

```{r example1}
discovery=GMHI[[1]];
validation=GMHI[[2]];
res=MRS(discovery, validation, GroupID="Group", DA.method="ancombc", measurement="shannon")
AUC=res[[3]]
```
1.2 using ALDEx2 method and Shannon index

```{r example2}
res=MRS(discovery, validation, GroupID="Group", DA.method="ALDEx2", measurement="shannon")
AUC=res[[3]]
```

### Example 2 Evaluation of MRS in terms of comparison between Healthy and a specific disease Healthy vs. CA

Calculate p-values based on 200 times of permutations in a single random data split 

```{r example3}
discovery.sub=prune_samples(sample_data(discovery)$Group1 %in% c("Healthy","CA"),discovery)
validation.sub=prune_samples(sample_data(validation)$Group1 %in% c("Healthy","CA"),validation)
res=MRS(discovery.sub, validation.sub, GroupID="Group", DA.method="ALDEx2", measurement="shannon")
AUC=res[[3]]
```


Development team
-------------
### Authors/developers
*MRS* is developed by:

* Chan Wang
* [Huilin Li](https://sites.google.com/site/huilinli09/)


### Maintainer
* Chan Wang 
* [Huilin Li](https://sites.google.com/site/huilinli09/)
