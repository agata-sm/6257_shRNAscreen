#! /usr/bin/env nextflow


/* 
 * minimal pipeline consisting of fastqc and multiqc
 * Using Nextflow DSL2
 * 
 * Author: Agata Smialowska
 * October 2022
 */ 

nextflow.enable.dsl=2

params.pipelinename="Minimal QC pipeline"

/* 
 * pipeline input parameters 
 */
params.resdir = "results"
params.projdir = "$launchDir/${params.projname}"

params.outdir = "${params.projdir}/${params.resdir}"

params.logdir = 'logs'

//params.fastqR1="$params.fastqdir/*R1*fastq.gz"

params.fastqPE="$params.fastqdir/*_R{1,2}_001.fastq.gz"
params.fastq="$params.fastqdir/*fastq.gz"



log.info """\
 M I N I M A L   Q C - N F   P I P E L I N E
 ===================================
 
 fastq files directory: ${params.fastqdir}

 data         : ${params.fastqdir}
 outdir       : ${params.outdir}
 """
 .stripIndent()

println ""

/////////////////////////////
// process metadata files

// get the files and sample names
// println "Samples"
// filesf = file("$params.sampleinfo")
// filesf.withReader {
//     String line
//     while( line = it.readLine() ) {
//         println line

//     }
// }


// // get the list of contrasts and samples for mageck
// println "Comparisons"
// comparisonf = file("$params.comparisons")
// comparisonf.withReader {
//     String line
//     while( line = it.readLine() ) {
//         println line

//     }
// }

println ""
println ""

/////////////////////////////
// channels

// samples channel
// smpls_ch= Channel.fromPath(params.sampleinfo, checkIfExists:true)
// 	smpls_ch
// 	    .splitCsv(header:true, sep: '\t', strip: true)
// 	    .map{ (it.sample) }
// 	    .toList()
// 	    .toListString().replace(/[/,"").replace(/]/,"").replace(/ /,"")
// 	    //.view()
// 	    .set { smpls_ch }


// fastq file paths channel - list of paths
// fastqr1_ch= Channel.fromPath(params.sampleinfo, checkIfExists:true)
// 	fastqr1_ch
// 	    .splitCsv(header:true, sep: '\t', strip: true)
// 	    .map{ (it.file) }
// 	    .collect { "${params.fastqdir}/$it" }
// 	    //.view()
// 	    .set { fastqr1_ch }

// fastq file paths pairs channel - paths
//fastqr1_chPE= Channel.fromPath(params.fastqPE , checkIfExists:true)
//	fastqr1_chPE
	    //.view()
//	    .set { fastqr1_chPE }

read_pairs = Channel.fromFilePairs(params.fastqPE, checkIfExists: true )
	read_pairs
	    .view()
	    .set { read_pairs }


// fastq file paths channel - paths
fastqr1_ch2= Channel.fromPath(params.fastq , checkIfExists:true)
	fastqr1_ch2
	    //.view()
	    .set { fastqr1_ch2 }


// fa for index
fa_ch=Channel.fromPath(params.shLibraryFa , checkIfExists:true)

/////////////////////////////
// processes
include { fastqc; multiqc; idx; mapPE } from './shRNAproc-modules.nf'



/////////////////////////////
// workflows

//default reads

workflow {

	//index
	idx(fa_ch)

	//read processing
	idx_bowtie_ch=idx.out.idx_bowtie_ch
		idx_bowtie_ch
			.flatten()
			.collect()
			.view()
			.set{ idx_bowtie_ch }


	map_readsPE_ch=read_pairs
		map_readsPE_ch
			.combine(idx_bowtie_ch)
			//.view()
			.set {map_readsPE_ch}

//[M47, [/proj/snic2022-23-410/private/data/M47_R1_001.fastq.gz, /proj/snic2022-23-410/private/data/M47_R2_001.fastq.gz], /crex/proj/snic2022-23-410/nobackup/private/nbis6257/analysis/tst2/work/cc/a76d2c137ebba7b4661f2df6ce6d78/shRNA_Idx_bowtie2.1.bt2, /crex/proj/snic2022-23-410/nobackup/private/nbis6257/analysis/tst2/work/cc/a76d2c137ebba7b4661f2df6ce6d78/shRNA_Idx_bowtie2.2.bt2, /crex/proj/snic2022-23-410/nobackup/private/nbis6257/analysis/tst2/work/cc/a76d2c137ebba7b4661f2df6ce6d78/shRNA_Idx_bowtie2.3.bt2, /crex/proj/snic2022-23-410/nobackup/private/nbis6257/analysis/tst2/work/cc/a76d2c137ebba7b4661f2df6ce6d78/shRNA_Idx_bowtie2.4.bt2, /crex/proj/snic2022-23-410/nobackup/private/nbis6257/analysis/tst2/work/cc/a76d2c137ebba7b4661f2df6ce6d78/shRNA_Idx_bowtie2.rev.1.bt2, /crex/proj/snic2022-23-410/nobackup/private/nbis6257/analysis/tst2/work/cc/a76d2c137ebba7b4661f2df6ce6d78/shRNA_Idx_bowtie2.rev.2.bt2]


	mapPE(map_readsPE_ch)
	
	//mapPE(idx_bowtie_ch, read_pairs)
	//mapPE(read_pairs)

	//QC
	fastqc(fastqr1_ch2)

	//multiQC
	multiqc_ch=fastqc.out.fastqc_report_ch
	multiqc(multiqc_ch.collect())

}


