// modules for shRNA-proc.nf


// outdirs

params.fastqc="FastQC"
params.fastqcOut="${params.outdir}/${params.fastqc}"

params.multiqc="MultiQC"
params.multiqcOut="${params.outdir}/${params.multiqc}"

params.idx="Bowtie2-idx"
params.idxOut="${params.outdir}/${params.idx}"

params.map="mappedPE"
params.mapOut="${params.outdir}/${params.map}"







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
    tuple path('shRNA_Idx_bowtie2*bt2'), path('shRNA_Idx_bowtie2.rev*bt2'), emit: idx_bowtie_ch

    script:
    """
    bowtie2-build -f $shFasta shRNA_Idx_bowtie2
    """
}


process mapPE {
    publishDir params.mapOut, mode:'copy'
    label 'mid_mem'

    input:
    tuple val(pair_id), path(reads), path(idx_bowtie_ch)

    output:
    path "${pair_id}.mapped.bowtie2.bam"

    script:
    """
    bowtie2 -p 4 -a --very-sensitive --dovetail --fr -x $idx_bowtie_ch -q -1 $reads[0] -2 $reads[1]  | samtools view -hbo ${pair_id}.mapped.bowtie2.bam -
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







