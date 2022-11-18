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
	    //.view()
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
include { fastqc; multiqc; idx; trim_readsPE; mapPE } from './shRNAproc-modules.nf'



/////////////////////////////
// workflows

//default reads

workflow {

	//index
	idx(fa_ch)

	//read processing
	idx_bowtie_ch=idx.out.idx_bowtie_ch
		idx_bowtie_ch
			//.flatten()
			.collect()
			.map{[it]}
			//.view()
			.set{ idx_bowtie_ch }


	trim_readsPE(read_pairs)

	map_readsPE_ch=trim_readsPE.out.trimmed_reads_PE_ch
		map_readsPE_ch
			//.map{[it]}
			//.view()
			.combine(idx_bowtie_ch)
			.view()
			//.flatten()
			//.view()
			.set {map_readsPE_ch}


//[OVCAR3-ct1, /crex/proj/snic2022-23-410/nobackup/private/nbis6257/analysis/tst2/work/a3/528c6624b366de33f193bdce07b0b4/OVCAR3-ct1.trimmed_merged_R1.fastq, /crex/proj/snic2022-23-410/nobackup/private/nbis6257/analysis/tst2/work/a3/528c6624b366de33f193bdce07b0b4/OVCAR3-ct1.trimmed_merged_R2.fastq]]

	mapPE(map_readsPE_ch)
	
	//mapPE(idx_bowtie_ch, read_pairs)
	//mapPE(read_pairs)

	//QC
	fastqc(fastqr1_ch2)

	//multiQC
	multiqc_ch=fastqc.out.fastqc_report_ch
	multiqc(multiqc_ch.collect())

}


