# TENxXeniumData
The `TENxXeniumData` package aims to contain a collection of Xenium spatial 
transcriptomics datasets provided by 10x Genomics. These datasets have been 
formatted into the Bioconductor classes, the SpatialExperiment or 
SpatialFeatureExperiment (SFE), to facilitate seamless integration into 
various applications, including examples, demonstrations, and tutorials. 
The constructed data objects include gene expression profiles, per-transcript 
location data, centroid, segmentation boundaries 
(e.g., cell or nucleus boundaries), and image. More datasets will be added 
in later versions.

# Installation 
To install this package, start R (version "4.4") and enter:
```{r}
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

# The following initializes usage of Bioc devel
BiocManager::install(version='devel')

BiocManager::install("TENxXeniumData")

```
