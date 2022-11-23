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

params.cnttab="count_table"
params.cnttabOut="${params.outdir}/${params.cnttab}"

params.mageck="mageck_rra"
params.mageckOut="${params.outdir}/${params.mageck}"

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
    tuple val(pair_id), path("${pair_id}.trimmed_merged_R1.fastq"), path("${pair_id}.trimmed_merged_R2.fastq"), emit: trimmed_reads_PE_ch
    path("${pair_id}.cutadapt_trim_fwd.log")
    path("${pair_id}.cutadapt_trim_rc.log")

    script:
    """
    
    cutadapt -e 0.1 -O 30 -m 48 -M 50  -a $params.adastring_fwd -A $params.adastring_rc --pair-filter=both \
    --untrimmed-output ${pair_id}.noada.trimFwd.r1.fastq --untrimmed-paired-output ${pair_id}.noada.trimFwd.r2.fastq \
    --too-long-output ${pair_id}.toolong.trimFwd.r1.fastq --too-long-paired-output ${pair_id}.toolong.trimFwd.r2.fastq \
    -o ${pair_id}.trimFwd.r1.fastq -p ${pair_id}.trimFwd.r2.fastq $reads >${pair_id}.cutadapt_trim_fwd.log 2>&1

    cat ${pair_id}.noada.trimFwd.r1.fastq ${pair_id}.toolong.trimFwd.r1.fastq >${pair_id}.untrimmed.trimFwd.r1.fastq
    cat ${pair_id}.noada.trimFwd.r2.fastq ${pair_id}.toolong.trimFwd.r2.fastq >${pair_id}.untrimmed.trimFwd.r2.fastq

    cutadapt -e 0.1 -O 30 -m 48 -M 50  -A $params.adastring_fwd -a $params.adastring_rc --pair-filter=both \
    --untrimmed-output ${pair_id}.noada.trimRc.r1.fastq --untrimmed-paired-output ${pair_id}.noada.trimRc.r2.fastq \
    --too-long-output ${pair_id}.toolong.trimRc.r1.fastq --too-long-paired-output ${pair_id}.toolong.trimRc.r2.fastq \
    -o ${pair_id}.trimRc.r1.fastq -p ${pair_id}.trimRc.r2.fastq ${pair_id}.untrimmed.trimFwd.r1.fastq ${pair_id}.untrimmed.trimFwd.r2.fastq >${pair_id}.cutadapt_trim_rc.log 2>&1

    cat ${pair_id}.trimFwd.r1.fastq ${pair_id}.trimRc.r1.fastq > ${pair_id}.trimmed_merged_R1.fastq
    cat ${pair_id}.trimFwd.r2.fastq ${pair_id}.trimRc.r2.fastq > ${pair_id}.trimmed_merged_R2.fastq

    """
}


process mapPE {
    publishDir params.mapOut, mode:'copy'
    label 'mid_mem'

    input:
    tuple val(pair_id), path(r1), path(r2), path(idx_bowtie_ch)

    output:
    tuple val(pair_id), path("${pair_id}.mapped.bowtie2.bam"), emit: mappedPE_ch

    script:
    """
    bowtie2 -p 4 -a --very-sensitive --dovetail --fr -x shRNA_Idx_bowtie2 -q -1 $r1 -2 $r2  | samtools view -hbo ${pair_id}.mapped.bowtie2.bam -
    """
}

process filter_reads {
    publishDir params.filtOut, mode:'copy'
    label 'small'

    input:
    tuple val(pair_id), path(bam_unfilt)

    output:
    path "${pair_id}.mapped.filt_mapq255_NM${params.nm}.bowtie2.bam", emit: filtered_ch
    path "${pair_id}.read_stats.log"


    script:
    """
    module load bioinfo-tools
    module load samtools/1.8
    module load NGSUtils/0.5.9

    echo "all alignments in sample ${pair_id}" >${pair_id}.read_stats.log
    samtools view -f 64 $bam_unfilt | wc -l >>${pair_id}.read_stats.log

    echo "aligned read pairs" >>${pair_id}.read_stats.log
    samtools view -F 256 -f 2 -q 255 -hbo ${pair_id}.mapped.filt_mapq255.bam $bam_unfilt
    
    echo "alignments with MAPQ 255" >>${pair_id}.read_stats.log

    samtools view -f 64 ${pair_id}.mapped.filt_mapq255.bam | wc -l  >>${pair_id}.read_stats.log

    echo "applying filter NM ${params.nm}"
    bamutils filter ${pair_id}.mapped.filt_mapq255.bam ${pair_id}.mapped.filt_mapq255_${params.nm}.bowtie2.bam -mismatch ${params.nm} -properpair >>${pair_id}.read_stats.log
    """
}


process count_table {
    publishDir params.cnttabOut, mode:'copy'
    label 'small'

    input:
    path(bam_filt)

    output:
    path "${params.projname}.counts", emit: count_table_ch
    path "${params.projname}.counts.summary"


    script:
    """
    featureCounts -p --countReadPairs --fracOverlap 0.958 -F SAF -a ${params.annot} -o ${params.projname}.counts $bam_filt
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







