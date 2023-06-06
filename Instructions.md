# Instructions

This file contains instructions on how to run the pipeline on Rackham (Uppmax). 
Running it on another HPC or a local system may require changes (please contact me for details).

## Pipeline

The pipeline is installed at `/proj/snic2022-23-410/nobackup/private/nbis6257/6257_shRNAscreen/shRNA-proc/`.

## Resources

Software used by the pipeline is containerised. All containers are built on Docker and used via Singularity / Apptainer. 
Their location on Rackham is:

`/proj/snic2022-23-410/nobackup/private/nbis6257/analysis/containers`

For some of the containers, the path to the containr was explicitely specified in `nextflow.config`; 
Before using the pipeline on another system please adapt `nextflow.config`, section `profile.singularity`.

## Input data

`Fastq` file names need to follow convention `SAMPLE_R{1/2}_001.fastq.gz`, 
where SAMPLE is sample ID which identifies given sample in the output files and the report.

The links to fastq files with name change when appropriate are at `/proj/snic2022-23-410/private/data`

## Usage

Project / run specific pipeline variables are set in file `shRNAproc.config`.
This file is present in each pipeline run directory.

Example contents:

```
	params {
	        allocation = "snic2022-22-783"

	        projname = "6257_shRNA_proc_OVCAR8"

	        fastqdir = "/proj/snic2022-23-410/private/data/OVCAR8"

	        shLibraryFa= "/proj/snic2022-23-410/nobackup/private/nbis6257/analysis/reference/shRNAlibrary/SHPH01_LPhs_1-10.desc.fa"

	        annot= "/proj/snic2022-23-410/nobackup/private/nbis6257/analysis/reference/shRNAlibrary/SHPH01_LPhs_1-10.saf"

	        // adapters
			    ada5="CGATTTCTTGGCTTTATATATCTTGTGGAAAGGACGAAACACCGG"
			    ada3="TTTTTGAATTCTCGACCTCGAGACAAATGGCAGTATTCATCCAC"

			    ada5rc="GTGGATGAATACTGCCATTTGTCTCGAGGTCGAGAATTCAAAAA"
			    ada3rc="CCGGTGTTTCGTCCTTTCCACAAGATATATAAAGCCAAGAAATCG"

		      // number of mismatches to tolerate in alignments
			    nm="1"
	}
```


We use program `screen` to isolate the pipeline session from the login session. 
This is required to allow the nextflow process to complete, i.e. the full analysis to be run till the end.

Basic instruction on how to use `screen` can be found in file `training-3iii2023.md`.

From *inside the new screen session*

```
	module load java/OracleJDK_11.0.9
	module load bioinfo-tools
	module load Nextflow/22.10.1

	export APPTAINERENV_TMPDIR="/proj/snic2022-23-410/nobackup/private/nbis6257/analysis/containers"
	export NXF_HOME="/proj/snic2022-23-410/nobackup/private/nbis6257/analysis/nf"
  
	nextflow run /proj/snic2022-23-410/nobackup/private/nbis6257/6257_shRNAscreen/shRNA-proc/shRNAproc.nf -profile cluster,singularity -c shRNAproc.config
```

Now `screen` can be disconnected, and the just started nextflow process will run in the background and continue to submit jobs to the queue in the cluster.
