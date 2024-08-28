/*
* WORKFLOW TO QUANTIFY VIRAL READS FROM BULK RNASEQ
*/

// PARAMETERS
params.fastq = null
params.outdir = null
params.strandedness = null
params.readType = null
params.viral_fasta = null
params.viral_t2g = null
params.cdna = null
params.genome = null
params.gtf = null


// CREATE VIRAL KALLISTO INDEX
process INDEX_VIRAL {
    container 'https://depot.galaxyproject.org/singularity/kb-python:0.28.2--pyhdfd78af_2'
    cpus 4
    memory '16 GB'

    input:
    path cdna
    path genome
    path viralT2G
    path viralFASTA

    output:
    path 'viral.idx'

    script:
    """
    cat ${cdna} ${genome} > combined_cdna_dna.fa.gz

    kb ref \
    --aa --d-list combined_cdna_dna.fa.gz -i viral.idx \
    -g ${viralT2G} -t ${task.cpus} --workflow custom \
    ${viralFASTA}
    """
}

// CREATE HOST ORGANISM KALLISTO INDEX
process INDEX_HOST {
    container 'https://depot.galaxyproject.org/singularity/kb-python:0.28.2--pyhdfd78af_2'
    cpus 4
    memory '16 GB'

    input:
    path cdna
    path genome
    path gtf

    output:
    path 'host.idx', emit: index
    path 'host_t2g.txt', emit: t2g

    script:
    """
    kb ref \
    -t ${task.cpus} -i host.idx -g host_t2g.txt \
    -f1 ${cdna} ${genome} ${gtf}
    """
}

// QUANIFY VIRAL READS
process QUANT_VIRAL {
    cpus 8
    memory '16 GB'
    container 'https://depot.galaxyproject.org/singularity/kb-python:0.28.2--pyhdfd78af_2'
    tag "${sampleID}"
    publishDir "${params.outdir}/viral/${sampleID}", mode: "copy"
    
    input:
    tuple val(sampleID), path(fastq)
    path viralIndex
    path viralT2G
    val strandedness
    val readType

    output:
    path 'quant_unfiltered/*.tsv'
    path '*.json'

    script:
    """
    kb count \
    -x BULK \
    -i ${viralIndex} \
    -g ${viralT2G} \
    --aa \
    --parity ${readType} \
    --strand ${strandedness} \
    --tcc \
    --matrix-to-files \
    -t ${task.cpus} \
    ${fastq}
    """
}

// QUANTIFY HOST READS
process QUANT_HOST {
cpus 8
    memory '16 GB'
    container 'https://depot.galaxyproject.org/singularity/kb-python:0.28.2--pyhdfd78af_2'
    tag "${sampleID}"
    publishDir "${params.outdir}/host/${sampleID}", mode: "copy"
    
    input:
    tuple val(sampleID), path(fastq)
    path hostIndex
    path hostT2G
    val strandedness
    val readType

    output:
    path 'quant_unfiltered/abundance*.tsv'
    path '*.json'

    script:
    """
    kb count \
    -x BULK \
    -i ${hostIndex} \
    -g ${hostT2G} \
    --parity ${readType} \
    --strand ${strandedness} \
    --tcc \
    --matrix-to-files \
    -t ${task.cpus} \
    ${fastq}
    """
}

// WORKFLOW
workflow {
    // import reads
    reads_ch = Channel.fromFilePairs(params.fastq, checkIfExists: true)
    // indexing
    INDEX_VIRAL(params.cdna, params.genome, params.viral_t2g, params.viral_fasta)
    INDEX_HOST(params.cdna, params.genome, params.gtf)
    // quantification
    QUANT_VIRAL(reads_ch, INDEX_VIRAL.out, params.viral_t2g, params.strandedness, params.readType)
    QUANT_HOST(reads_ch, INDEX_HOST.out.index, INDEX_HOST.out.t2g, params.strandedness, params.readType)
}