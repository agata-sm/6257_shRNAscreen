// modules for shRNA-proc.nf


// outdirs

params.fastqc="FastQC"
params.fastqcOut="${params.outdir}/${params.fastqc}"

params.multiqc="MultiQC"
params.multiqcOut="${params.outdir}/${params.multiqc}"

params.trim="trimmedPE"
params.trimOut="${params.outdir}/${params.trim}"


params.idx="Bowtie2-idx"
params.idxOut="${params.outdir}/${params.idx}"

params.map="mappedPE"
params.mapOut="${params.outdir}/${params.map}"

params.filt="filteredPE"
params.filtOut="${params.outdir}/${params.filt}"

params.filtLogs="read_logs"
params.filtLogsOut="${params.outdir}/${params.filtLogs}"

params.cnttab="count_table"
params.cnttabOut="${params.outdir}/${params.cnttab}"

params.cnttabProc="count_table_processed"
params.cnttabProcOut="${params.outdir}/${params.cnttabProc}"

// format adapter strings
params.adastring_fwd="${params.ada5}...${params.ada3}"
params.adastring_rc="${params.ada5rc}...${params.ada3rc}"



// scripts
params.scripts="${projectDir}/bin"

// versions
params.verfile="software.versions"

/// modules


process fastqc {
    publishDir params.fastqcOut, mode:'copy'
    label 'small'

    input:
    path fastqfile

    output:
    path('*'), emit: fastqc_report_ch
    //path "${params.verfile}"

    script:
    """
    echo "fastqc $fastqfile"
    fastqc $fastqfile

    #echo "Software versions for ${params.pipelinename}" >${params.verfile}
    #date >>${params.verfile}
    #echo "process ** fastqc **" >>${params.verfile}
    #fastqc -v >>${params.verfile}
    """

}


process idx {
    publishDir params.idxOut, mode:'copy'
    label 'small'


    input:
    path shFasta

    output:
    //tuple path('shRNA_Idx_bowtie2.1.bt2'),path('shRNA_Idx_bowtie2.2.bt2'),path('shRNA_Idx_bowtie2.3.bt2'),path('shRNA_Idx_bowtie2.4.bt2'),path('shRNA_Idx_bowtie2.rev.1.bt2'),path('shRNA_Idx_bowtie2.rev.2.bt2') , emit: idx_bowtie_ch
    tuple path('shRNA_Idx_bowtie2.{1,2,3,4}.bt2'),path('shRNA_Idx_bowtie2.rev.{1,2}.bt2') , emit: idx_bowtie_ch

    script:
    """
    bowtie2-build -f $shFasta shRNA_Idx_bowtie2
    """
}


process trim_readsPE {
    publishDir params.trimOut, mode:'copy'
    label 'small'

    input:
    tuple val(pair_id), path(reads)

    output:
    tuple val(pair_id), path("${pair_id}.trimmed_merged_R1.fastq.gz"), path("${pair_id}.trimmed_merged_R2.fastq.gz"), emit: trimmed_reads_PE_ch
    path("${pair_id}.cutadapt_trim_fwd.log"), emit: trimlog_f_ch
    path("${pair_id}.cutadapt_trim_rc.log"), emit: trimlog_r_ch

    script:
    """
    
    cutadapt -e 0.1 -O 30 -m 48 -M 50  -a $params.adastring_fwd -A $params.adastring_rc --pair-filter=both \
    --untrimmed-output ${pair_id}.noada.trimFwd.r1.fastq.gz --untrimmed-paired-output ${pair_id}.noada.trimFwd.r2.fastq.gz \
    --too-long-output ${pair_id}.toolong.trimFwd.r1.fastq.gz --too-long-paired-output ${pair_id}.toolong.trimFwd.r2.fastq.gz \
    -o ${pair_id}.trimFwd.r1.fastq.gz -p ${pair_id}.trimFwd.r2.fastq.gz $reads >${pair_id}.cutadapt_trim_fwd.log 2>&1

    cat ${pair_id}.noada.trimFwd.r1.fastq.gz ${pair_id}.toolong.trimFwd.r1.fastq.gz >${pair_id}.untrimmed.trimFwd.r1.fastq.gz
    cat ${pair_id}.noada.trimFwd.r2.fastq.gz ${pair_id}.toolong.trimFwd.r2.fastq.gz >${pair_id}.untrimmed.trimFwd.r2.fastq.gz

    cutadapt -e 0.1 -O 30 -m 48 -M 50  -A $params.adastring_fwd -a $params.adastring_rc --pair-filter=both \
    --untrimmed-output ${pair_id}.noada.trimRc.r1.fastq.gz --untrimmed-paired-output ${pair_id}.noada.trimRc.r2.fastq.gz \
    --too-long-output ${pair_id}.toolong.trimRc.r1.fastq.gz --too-long-paired-output ${pair_id}.toolong.trimRc.r2.fastq.gz \
    -o ${pair_id}.trimRc.r1.fastq.gz -p ${pair_id}.trimRc.r2.fastq.gz \
    ${pair_id}.untrimmed.trimFwd.r1.fastq.gz ${pair_id}.untrimmed.trimFwd.r2.fastq.gz >${pair_id}.cutadapt_trim_rc.log 2>&1

    cat ${pair_id}.trimFwd.r1.fastq.gz ${pair_id}.trimRc.r1.fastq.gz > ${pair_id}.trimmed_merged_R1.fastq.gz
    cat ${pair_id}.trimFwd.r2.fastq.gz ${pair_id}.trimRc.r2.fastq.gz > ${pair_id}.trimmed_merged_R2.fastq.gz

    """
}


process mapPE {
    publishDir params.mapOut, mode:'copy'
    label 'big_mem'
    cpus params.threads_bigmem
    scratch true


    input:
    tuple val(pair_id), path(r1), path(r2), path(idx_bowtie_ch)

    output:
    tuple val(pair_id), path("${pair_id}.mapped.bowtie2.bam"), emit: mappedPE_ch

    script:
    """
    bowtie2 -p ${params.threads_bigmem} --quiet -a --very-sensitive --dovetail --fr -x shRNA_Idx_bowtie2 -q -1 $r1 -2 $r2  | samtools view -hbo ${pair_id}.mapped.bowtie2.bam -
    """
}

process filter_alns {
    publishDir params.filtOut, mode:'copy'
    label 'small'

    input:
    tuple val(pair_id), path(bam_unfilt)

    output:
    path "${pair_id}.mapped.filt_mapq255_NM${params.nm}.bowtie2.bam", emit: filtered_ch
    path "${pair_id}.read_stats.log", emit: readlogs_ch


    script:
    """
    #module load bioinfo-tools
    #module load samtools/1.8
    #module load NGSUtils/0.5.9

    echo "all R1 alignments (proxy for aligned read pairs) in sample ${pair_id}" >${pair_id}.read_stats.log
    echo "aligned read pairs: `samtools view -f 64 $bam_unfilt | wc -l`" >>${pair_id}.read_stats.log
    samtools view -f 64 $bam_unfilt | wc -l >>${pair_id}.read_stats.log
    echo "" >>${pair_id}.read_stats.log

    samtools view -F 256 -f 2 -q 255 -hbo ${pair_id}.mapped.filt_mapq255.bam $bam_unfilt
    
    echo "R1 alignments with MAPQ 255 (proxy for correctly aligned read pairs)" >>${pair_id}.read_stats.log
    echo "aligned read pairs with MAPQ 255: `samtools view -f 64 ${pair_id}.mapped.filt_mapq255.bam | wc -l `" >>${pair_id}.read_stats.log
    samtools view -f 64 ${pair_id}.mapped.filt_mapq255.bam | wc -l  >>${pair_id}.read_stats.log

    echo "" >>${pair_id}.read_stats.log
    echo "applying filter NM ${params.nm}; reads counted ">>${pair_id}.read_stats.log
    bamutils filter ${pair_id}.mapped.filt_mapq255.bam ${pair_id}.mapped.filt_mapq255_NM${params.nm}.bowtie2.bam -mismatch ${params.nm} -properpair >>${pair_id}.read_stats.log
    
    echo "" >>${pair_id}.read_stats.log
    echo "R1 alignments with MAPQ 255 which passed the mismatch filter (proxy for aligned read pairs with max ${params.nm} mismatches)" >>${pair_id}.read_stats.log
    echo "aligned read pairs with MAPQ 255 filtered for NM ${params.nm}: `samtools view -f 64 ${pair_id}.mapped.filt_mapq255_NM${params.nm}.bowtie2.bam | wc -l `">>${pair_id}.read_stats.log
    samtools view -f 64 ${pair_id}.mapped.filt_mapq255_NM${params.nm}.bowtie2.bam | wc -l >>${pair_id}.read_stats.log

    """
}


process read_logs {
    publishDir params.filtLogsOut, mode:'copy'
    label 'small'

    input:
    path(readlogs_ch)
    path(trimlog_f_ch)
    path(trimlog_r_ch)

    output:
    path "log_stats.txt"

    script:

    """
    perl ${params.scripts}/6257_parse_logs.pl --outdir .
    """

}

process count_table {
    publishDir params.cnttabOut, mode:'copy'
    label 'mid_mem'

    input:
    path(bam_filt)

    output:
    path "${params.projname}.counts", emit: count_table_ch
    path "${params.projname}.counts.summary"


    script:
    """
    #module load bioinfo-tools
    #module load subread/2.0.3
    
    featureCounts -p -B -C --fracOverlap 0.958 -F SAF -a ${params.annot} -o ${params.projname}.counts $bam_filt
    """
}

process filt_count_table {
    publishDir params.cnttabProcOut, mode:'copy'
    label 'small'

    input:
    path(count_table_ch)

    output:
    path "${params.projname}.counts_processed.all.tsv"
    path "${params.projname}.counts_processed.0rm.tsv"
    path "${params.projname}.counts_processed.0rm_noAlt.tsv"


    script:
    """
    #module load bioinfo-tools
    #module load perl/5.26.2

    perl ${params.scripts}/6257_proc_cnt_table.pl --infile ${params.projname}.counts --outfile  ${params.projname}.counts_processed.all.tsv --library ${params.libraryDescription} --setup all

    perl ${params.scripts}/6257_proc_cnt_table.pl --infile ${params.projname}.counts --outfile  ${params.projname}.counts_processed.0rm.tsv --library ${params.libraryDescription} --setup 0rm

    perl ${params.scripts}/6257_proc_cnt_table.pl --infile ${params.projname}.counts --outfile  ${params.projname}.counts_processed.0rm_noAlt.tsv --library ${params.libraryDescription} --setup 0rm_noAlt
    """
}



process multiqc {
    publishDir params.multiqcOut, mode:'copy'
    label 'small'


    input:
    path fastqc_report
    
    output:
    path "multiqc_report.html"
    path "multiqc_data"

    script:
    """
    multiqc .
    """
}







