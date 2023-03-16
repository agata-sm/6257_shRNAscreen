#! /usr/bin/env nextflow


/* 
 * minimal pipeline consisting of fastqc and multiqc
 * Using Nextflow DSL2
 * 
 * Author: Agata Smialowska
 * October 2022
 */ 

nextflow.enable.dsl=2

params.pipelinename="shRNA Screen Processing"

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
 SH-RNA SCREEN PROCESSING - N F   P I P E L I N E
 ==================================================
 
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
include { fastqc; multiqc; idx; trim_readsPE; mapPE; filter_alns; count_table; filt_count_table; read_logs } from './shRNAproc-modules.nf'



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
			.map{[it]}
			//.view()
			.set{ idx_bowtie_ch }


	trim_readsPE(read_pairs)

	map_readsPE_ch=trim_readsPE.out.trimmed_reads_PE_ch
		map_readsPE_ch
			.combine(idx_bowtie_ch)
			//.view()
			.set {map_readsPE_ch}


	mapPE(map_readsPE_ch)
	
	filter_alns(mapPE.out.mappedPE_ch)

	filt_bams_ch=filter_alns.out.filtered_ch
	count_table(filt_bams_ch.collect())

	//trimlog_f_ch=trim_readsPE.out.trimlog_f_ch.collect()
	//trimlog_r_ch=trim_readsPE.out.trimlog_r_ch.collect()
	filt_bams_logs_ch=filter_alns.out.readlogs_ch
	read_logs(filt_bams_logs_ch.collect(), trim_readsPE.out.trimlog_f_ch.collect(), trim_readsPE.out.trimlog_r_ch.collect())
	
	filt_count_table(count_table.out.count_table_ch)

	//QC
	fastqc(fastqr1_ch2)

	//multiQC
	multiqc_ch=fastqc.out.fastqc_report_ch
	multiqc(multiqc_ch.collect())

}


