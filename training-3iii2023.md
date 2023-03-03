# Pipeline usage

10 Mar 2023


## Software requirements

The pipeline for processing shRNA sequencing reads is designed to work on HPC clusters which provide tools for efficient process parallelisation (such as Rackham). It uses Docker / Singularity (Apptainer) containers to ensure results reproducibility. The output of the pipeline are summarised read pairs mapped to shRNA hairpin sequences in the form of count tables. 

Besides running the data processing on an HPC cluster, the final count tables need to be further analysed to compute summary statistics, detect enriched shRNAs, and finally produce a report. The last step of this analysis, report compilation can be performed locally using R. To do this, the following software needs to be installed.

* recent version of R (https://cran.r-project.org/)

* RStudio (https://posit.co/downloads/)

* R packages::

```
install.packages("knitr")
install.packages("bookdown")

install.packages("tidyverse")
install.packages("ggrepel")
install.packages("reshape2")
install.packages("magrittr")
install.packages("cowplot")
install.packages("ggVennDiagram")

install.packages("viridis")
install.packages("DT")
install.packages("htmltools")

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("MAGeCKFlute","edgeR","org.Hs.eg.db")
```


 **Note**
 Please check that the software is installed before the training session, as installing all of it at once will take some time.
 
 

:computer:

:office:

:globe_with_meridians:
