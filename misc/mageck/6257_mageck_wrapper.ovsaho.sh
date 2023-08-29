#!/bin/bash -l
#SBATCH -A snic2022-22-783
#SBATCH -p core
#SBATCH -n 4
#SBATCH -t 14:00:00
#SBATCH -J mageck_ovsh
#SBATCH -o mageck.ovsaho.29viii2023.out
#SBATCH -e mageck.ovsaho.29viii2023.out


module load bioinfo-tools
module load MAGeCK/0.5.9.4
module load R_packages/4.1.1
module load pandoc/2.17.1.1





### directory where the count tables are saved
### usually in results/count_table_processed after shRNAproc.nf run
### files suffixed with the following will be analysed (standard pipeline output)
### PROJ.counts_processed.0rm_noAlt.tsv
### PROJ.counts_processed.0rm.tsv
### PROJ.counts_processed.all.tsv


indir="/proj/snic2022-23-410/nobackup/private/nbis6257/analysis/OVSAHO/6257_shRNA_proc_OVSAHO/results/count_table_processed/"


### identical to projname in nextflow run.config; PROJ in the file names above
project_prefix="6257_shRNA_proc_OVSAHO"



### outdir - location to save the results
outdir="/crex/proj/snic2022-23-410/nobackup/private/nbis6257/analysis/OVSAHO/mageck"


###### contrasts and samples

### assuming three comparisons:
### treated vs untreated
### treated vs ctrl
### untreated vs ctrl

### sample names have to be identical to header in counts table
### which are SMPL derived from fastq file names SMPL_R1/2_001.fastq

## contrast1

# name
contrast1="M19_24_vs_M13_18"
# samples in each group
control1="M13,M14,M15,M16,M17,M18"
treated1="M19,M20,M21,M22,M23,M24"


### contrast2

# name
contrast2="M13_18_vs_ctrl"
# samples in each group
control2="OVSAHOCT2,OVSAHOCT3"
treated2="M13,M14,M15,M16,M17,M18"


### contrast3

# name
contrast3="M19_24_vs_ctrl"
# samples in each group
control3="OVSAHOCT2,OVSAHOCT3"
treated3="M19,M20,M21,M22,M23,M24"





################################### end of input

mkdir -p $outdir


declare -A control_smpls
declare -A treatment_smpls

control_smpls[$contrast1]=$control1
control_smpls[$contrast2]=$control2
control_smpls[$contrast3]=$control3

treatment_smpls[$contrast1]=$treated1
treatment_smpls[$contrast2]=$treated2
treatment_smpls[$contrast3]=$treated3


for contrast in $contrast1 $contrast2 $contrast3
do
	echo "processing comparison $contrast"

	conrast_pref="${project_prefix}.${contrast}"

	outdir_contrast="${outdir}/${conrast_pref}"

	table_suf1="0rm_noAlt"
	table_suf2="0rm"
	table_suf3="all"

	for table_suf in $table_suf1 $table_suf2 $table_suf3
	do

		conrast_table_pref="${project_prefix}.${contrast}.${table_suf}"

		outdir_contrast_suf="${outdir_contrast}/${table_suf}"
		mkdir -p $outdir_contrast_suf


		count_table="${indir}/${project_prefix}.counts_processed.${table_suf}.tsv"

		cd $outdir_contrast_suf

		if [ $table_suf == "all" ]
		then
			echo "mageck test -k $count_table -c ${control_smpls[$contrast]} -t ${treatment_smpls[$contrast]} -n $conrast_table_pref --norm-method total --pdf-report"
			mageck test -k $count_table -c ${control_smpls[$contrast]} -t ${treatment_smpls[$contrast]} -n $conrast_table_pref --norm-method total --pdf-report
			echo ""

		else
			echo "mageck test -k $count_table -c ${control_smpls[$contrast]} -t ${treatment_smpls[$contrast]} -n $conrast_table_pref --norm-method total --pdf-report --remove-zero both"
			mageck test -k $count_table -c ${control_smpls[$contrast]} -t ${treatment_smpls[$contrast]} -n $conrast_table_pref --norm-method total --pdf-report --remove-zero both
			echo ""

		fi
	done

	echo "** done comparison $contrast"
	echo ""

done

echo "*** ALL DONE ***"



