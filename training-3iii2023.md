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
 
 
 ## Data processing pipeline
 
This pipeline uses a workflow manager called `nextflow`. It allows for automated task parallelisation and management, while ensuring process reproducibility via Docker containers.
 
 More about `nextflow` can be found at [Nextflow documentation](https://www.nextflow.io/docs/latest/basic.html).
 
 More about Docker can be found at [Docker documentation](https://docs.docker.com/).

The structure of the shRNA processing pipeline is shown on Figure 1.


In practice, the workflow is started by a single command, and runs without user intervention.

### Before the run

The pipeline requires certain project specific variables set before the run. This is done by modyfying a **config** file in **project directory** where the pipeline is being started.

Example of a config file `shRNAproc.config` is given below.


```
	params {
	        allocation = "snic2022-22-634"

	        projname = "6257_shRNA_proc_test3"

	        fastqdir = "/proj/snic2022-23-333/nobackup/private/agata/nbis6257/analysis/test/data"

	        shLibraryFa= "/proj/snic2022-23-333/nobackup/private/agata/nbis6257/analysis/library/SHPH01_LPhs_1-10.desc.fa"

	        annot= "/proj/snic2022-23-333/nobackup/private/agata/nbis6257/analysis/library/SHPH01_LPhs_1-10.saf"

	        libraryDescription="/proj/snic2022-23-333/nobackup/private/agata/nbis6257/analysis/library/SHPH01_LPhs_1-10.desc2.nf.tsv"

	        // adapters
			ada5="CGATTTCTTGGCTTTATATATCTTGTGGAAAGGACGAAACACCGG"
			ada3="TTTTTGAATTCTCGACCTCGAGACAAATGGCAGTATTCATCCAC"

			ada5rc="GTGGATGAATACTGCCATTTGTCTCGAGGTCGAGAATTCAAAAA"
			ada3rc="CCGGTGTTTCGTCCTTTCCACAAGATATATAAAGCCAAGAAATCG"

	        // number of mismatches to tolerate in alignments
			nm="1"

	}
```


The contents:

* `allocation` is the computing allocation on Rackham - required to submit jobs to the cluster queue;
* `projname` is the project specific prefix added to result files and directories;
* `fastqdir` is a directory where paired - end **fastq** files are saved; please look below for file naming;
* `shLibraryFa` is fasta file with shRNA hp sequences;
* `annot` is annotation for read summarisation;
* `libraryDescription` is description of shRNA probes derived from library files.

* section **adapters** lists adapters used for library cloning;

* `nm` - number of mismatches to tolerate in alignments to pass the quality filter.

Example of this file can be found in directory `proj-config-files`; The locations of `shLibraryFa`, `annot` and `libraryDescription` are set for the current storage project on Rackham and can be used without change.

**Note**
Only `projname`  and its corresponding `fastqdir` need to be updated if `proj-config-files/shRNAproc.config` is used for pipeline run.

Please see NBIS Project Report for documentation on how `shLibraryFa`, `annot` and `libraryDescription` were obtained.

#### Comment on naming fastq files

Depending on the sequencing provider, naming of fastq files may follow different conventions. The processing pipeline requires that the files contain suffix `R1_001.fastq` (e.g. `M10_R1_001.fastq.gz`). Fastq files may be compressed (`.gz`) or not (gzipped files are recommended).


### Preparation

:office: :globe_with_meridians:
This part is run on a remote server, i.e. Rackham.
Login to Rackham first.

To start with, you need to change the directory to project directory:


```
cd /proj/snic2022-23-410/nobackup/private
```

Let's keep the analysis from other parts of the project:

```
cd analysis
```

You will create a directory for each pipeline run, and change to it, for example:


```
mkdir testrun1
cd testrun1
```



:computer:

:office:

:globe_with_meridians:
