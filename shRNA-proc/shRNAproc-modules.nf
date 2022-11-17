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
    path "${params.verfile}"

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
    path 'shRNA_Idx_bowtie2*', emit: idx_bowtie

    script:
    """
    bowtie2-build -f $shFasta shRNA_Idx_bowtie2
    """
}


process mapPE {
    publishDir params.mapOut, mode:'copy'
    label 'small'

    input:
    path idx_bowtie
    tuple val(pair_id), path(r1fq,r2fq)

    output:


    script:
    """
    echo $idx_bowtie
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







