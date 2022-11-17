// modules for minimal.nf


params.fastqc="FastQC"
params.fastqcOut="${params.outdir}/${params.fastqc}"

params.multiqc="MultiQC"
params.multiqcOut="${params.outdir}/${params.multiqc}"



// assets
//params.countertemplate="${projectDir}/assets/template.properties"


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
    fastqc $fastqr1

    echo "Software versions for ${params.pipelinename}" >${params.verfile}
    date >>${params.verfile}
    echo "process ** fastqc **" >>${params.verfile}
    fastqc -v >>${params.verfile}
    """

}


process multiqc {
    publishDir params.multiqcOut, mode:'copy'
    label 'small'


    input:
    file ('fastqc/*') from fastqc_report_ch

    output:
    file "multiqc_report.html" into multiqc_report
    file "multiqc_data"

    script:
    """
    multiqc .
    """
}






