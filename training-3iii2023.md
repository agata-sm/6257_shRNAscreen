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

* `git` (installation instructions at https://git-scm.com/downloads)


 **Note**
 Please check that the software is installed before the training session, as installing all of it at once will take some time.
 
 
 </br>
 </br>
 
 ## Data processing pipeline
 
This pipeline uses a workflow manager called `nextflow`. It allows for automated task parallelisation and management, while ensuring process reproducibility via Docker containers.
 
 More about `nextflow` can be found at [Nextflow documentation](https://www.nextflow.io/docs/latest/basic.html).
 
 More about Docker can be found at [Docker documentation](https://docs.docker.com/).

The structure of the shRNA processing pipeline is shown on Figure 1.


In practice, the workflow is started by a single command, and runs without user intervention.


[![](https://mermaid.ink/img/pako:eNp9VMtu4yAU_ZWIlSu1I4PtJM2imz52lSLN7OzKurFxjQQ4xUSdUdN_HzCGhjQzXvmec19wjvhAzdBStEEdH96bHpRe_Hqo5MJ8-zRJ7nuQkvIfnRrEE-N0C0yNV1dzAo4TtqD7wJH_cFnJ2t8vc5AnZcdBaypfPF8kZTNwThsdoGVSCtiHcFVqxUStKLTj9tG3WifJIgy5PQ1waluKHZM09MDYtvyqxqTsGNdU1cDZqxRU6jFwWdQt_74gLgx0kLrWsOM01C2jutU0ob6UGO2Oo-VJGkW47GDUb40vJeT7NiQrxYFrdpKVR02K0yhd3NzcHe1l1nsr8NFcsN9koqaBCtdNT452A99lJg1u4MyrO6FG4Xo3vGtGHZt7uS1rJPZSu3DpZb5Yi1Mv-0Qb1bzyl3gjfRyugwn-UY6DI1wBJsES7oRWtB2Icc7PI9q6IxhjBorgihlYnQPL4AkHkPQMwLfnwDo4wJ3DCex2IiS4Ye6XBSvMQHEOmGOgaySoEsBa8wR82IQK6Z4KWqGN-W1pB2ZKhSr5aVLhoIeff2SDNh3wkV6jw74FTR8YvCoQAaUt04N6du_K9Lx8_gVCbkIp?type=png)](https://mermaid.live/edit#pako:eNp9VMtu4yAU_ZWIlSu1I4PtJM2imz52lSLN7OzKurFxjQQ4xUSdUdN_HzCGhjQzXvmec19wjvhAzdBStEEdH96bHpRe_Hqo5MJ8-zRJ7nuQkvIfnRrEE-N0C0yNV1dzAo4TtqD7wJH_cFnJ2t8vc5AnZcdBaypfPF8kZTNwThsdoGVSCtiHcFVqxUStKLTj9tG3WifJIgy5PQ1waluKHZM09MDYtvyqxqTsGNdU1cDZqxRU6jFwWdQt_74gLgx0kLrWsOM01C2jutU0ob6UGO2Oo-VJGkW47GDUb40vJeT7NiQrxYFrdpKVR02K0yhd3NzcHe1l1nsr8NFcsN9koqaBCtdNT452A99lJg1u4MyrO6FG4Xo3vGtGHZt7uS1rJPZSu3DpZb5Yi1Mv-0Qb1bzyl3gjfRyugwn-UY6DI1wBJsES7oRWtB2Icc7PI9q6IxhjBorgihlYnQPL4AkHkPQMwLfnwDo4wJ3DCex2IiS4Ye6XBSvMQHEOmGOgaySoEsBa8wR82IQK6Z4KWqGN-W1pB2ZKhSr5aVLhoIeff2SDNh3wkV6jw74FTR8YvCoQAaUt04N6du_K9Lx8_gVCbkIp)

Figure 1. Structure of the shRNA-seq processing pipeline.

 </br>
 </br>

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

The locations of `shLibraryFa`, `annot` and `libraryDescription` are set for the current storage project on Rackham and can be used without change.

**Note**
Only `projname`  and its corresponding `fastqdir` need to be updated if `proj-config-files/shRNAproc.config` is used for pipeline run.

Please see NBIS Project Report for documentation on how `shLibraryFa`, `annot` and `libraryDescription` were obtained. Shortened examples are given in directory `misc/library`.


#### Comment on naming fastq files

Depending on the sequencing provider, naming of fastq files may follow different conventions. The processing pipeline requires that the files contain suffix `R1_001.fastq` and `R2_001.fastq` (e.g. `M10_R1_001.fastq.gz`). Fastq files may be compressed (`.gz`) or not (gzipped files are recommended though, to save space).

While this is not a strict requirement, I recommend to shorten file names to follow a convention `SMPL_R1/2_001.fastq.gz` where `SMPL` is sample name. This is because the part of file names prior to `_R1/2_001.fastq` will be used to name output files and as sample ID in the report. As a consequence, file names of style `P1234_FXNNB55578_R1_001.fastq` will result in hard to read reports. At the same time, it is considered best practice to keep the original file names of raw data. One way to address this seemingly contradictory situation without resorting to copying files under redacted names is to **create soft links** to fastq files.

For example, in desired location which will be used as `fastqdir` in the config file:

```
ln -s /path/to/data/P1234_FXNNB55578_R1_001.fq.gz SMPL1_R1_001.fastq.gz
```

Link `SMPL1_R1_001.fastq.gz` will point to the original file `path/to/data/P1234_FXNNB55578_R1_001.fq.gz`.

Please bear in mind that when creating links it is **essential** to create *soft links* by using `ln -s` and not just `ln` - the latter may lead to accidentally removing the original file and data loss.




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

To shorten the commands, you can create a variable which holds the path to the pipeline directory:

```
export PIPELINE_DIR="/proj/snic2022-23-333/nobackup/private/agata/nbis6257/tsts/6257_shRNAscreen/"
```


Copy the config file `shRNAproc.config` from the pipeline directory to your working directory:

```
cp $PIPELINE_DIR/shRNA-proc/shRNAproc.config .
```

You can now modify the config file to the needs of the current analysis. In this example we do this by using a simple text editor `nano`, but you can also prepare it beforehand on your local computer using text editor.

```
nano shRNAproc.config

### make edits 
### in this case change nothing as the path is set already to test data
### usually you need to chanhe the path to fastq files
```

To save file and quit, follow the instruction at the bottom of the screen. To quit is `Ctrl-X`, you get prompted if to save file, press `Y` and then `Enter` - we keep this current file name.

You can inspect the contents of this changed file, to confirm all is as it should be:


```
cat shRNAproc.config
```

Now we are almost ready to start the pipeline. 

### Running the shRNA seq processing pipeline

:office: :globe_with_meridians:
This part is run on a remote server, i.e. Rackham.
You are already logged to Rackham.

It is most practical to run the pipeline in the background. (In fact this is the only viable way.) This protects the run from accidental session interruption - for example when you connect remotely to the server and the session disconnects.

You can use several programs to achieve this, in this example we use `screen`, which is usually already installed in any Linux distribtion.

First, start the program by typing (on login node):

```
screen 
```

A new terminal appears. You can start a process in it, disconnect from it, then reconnect at any time (`screen -r`).

To start a new screen press `Ctrl-a`, then `c`. To run the pipeline:

```
module load java/OracleJDK_11.0.9
module load bioinfo-tools
module load Nextflow/22.10.1

APPTAINERENV_TMPDIR="/proj/snic2022-23-333/nobackup/private/agata/nbis6257/containers"

NXF_HOME="/proj/snic2022-23-333/nobackup/private/agata/nbis6257/analysis/nf"

nextflow run $PIPELINE_DIR/shRNA-proc/shRNAproc.nf -profile cluster,singularity -c shRNAproc.config
```

You will see messages printed on the screen, as the run progresses.

You can now disconnect from this process, by pressing `Ctrl-a` then `d`.

You can now monitor the progression by checking jobs submitted to the queue:

```
jobinfo -u $USER
```

The run on the complete data set will take several hours. 

 </br>
 </br>


## Statistical analysis using MAGeCK

:office: :globe_with_meridians:
This part is run on a remote server, i.e. Rackham.
Login to Rackham first and inspect whether the pipeline run has finished.

 </br>
 </br>


After the processing run will have completed, we can use its final output to perform statistical analysis. We will use [MAGeCK](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-014-0554-4) to perform robust rank aggregatoin (RRA) analysis. Detailed description of the program and its many options can be found at [MAGeCK homepage](/https://sourceforge.net/p/mageck/wiki/Home/]).

Script `6257_mageck_wrapper.sh` handles the software modules and the commands. You can copy it from `$PIPELINE_DIR`:

```
cp $PIPELINE_DIR/misc/mageck/6257_mageck_wrapper.sh .
```

 </br>
 </br>


For a new analysis you will need to change some details (directory locations, sample names, contrasts, computing project allocation), the file is ready for the test run.

Please bear in mind that **sample names have to be identical to samples in the header** of the counts tables `PROJ_PREFIX/results/count_table_processed/PROJ_PREFIX.counts_processed.*.tsv` where `PROJ_PREFIX` is `projname` in the config file. In practice, these names are `SMPL` derived from fastq file names `SMPL_R1/2_001.fastq`.

 </br>
 </br>


After modyfying the script you can submit it to the queue:

```
sbatch 6257_mageck_wrapper.sh
```

 </br>
 </br>

Inspect the queue as before 

```
jobinfo -u $USER
```

and after the run is over, you can also check the log file saved in the directory where the script was run. 

The results are ready for report compilation.

 </br>
 </br>


## Report

:computer:
This part is run on a local computer.

 </br>
 </br>


First, we need to copy the report data to the local computer.

We will need:

* Results of RRA from MAGeCK;

* processed count tables from `PROJ_PREFIX/results/count_table_processed/`

* metadata files describing samples `metadata.txt` and comparisons `comparisons.txt` (examples can be found at `misc/metadata`)

* read summarisation statistics at `PROJ_PREFIX/results/read_logs/log_stats.txt`

 </br>
 </br>


:office: :globe_with_meridians:
This part is run on a remote server, i.e. Rackham.

Login to Rackham first and inspect whether the pipeline run has finished. Copy the required directories using the method of choice.

 </br>
 </br>


:computer:
This part is run on a local computer.

You will need the report template, script and other files for report compilation. You can get them by cloning this repository and copying the files in `misc/report` to your working directory.


```
git clone https://github.com/agata-sm/6257_shRNAscreen.git
```

