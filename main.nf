#!/usr/bin/env nextflow

def helpMessage() {
    log.info """
    Usage:
    The typical command for running the pipeline is as follows:
    nextflow run main.nf --reference sample.fa --gtf sample.gtf [Options]
    
    Inputs Options:
    --reference         Reference file (fasta, gzipped or not)
    --gtf               Annotation file (GTF)

    Read Length Options:
      Use either:
    --read_length       Comma-separated list with read lengths (string)
                        (default: $params.read_length)
      Or:
    --read_length_from  Start value for read length (int)
                        (default: $params.read_length_from)
    --read_length_to    End value for read length (int)
                        (default: $params.read_length_to)
    --read_length_by    Value to increment the read length by (int)
                        If not provided, it defaults to an increment of 10.
                        (default: $params.read_length_by)
    
    STAR Options:
    --star_container    Nextflow-compatible container with STAR (string)
                        (default: $params.star_container)
    --index_n_bases     Value for STAR genomeSAindexNbases parameter (int)
                        (default: $params.index_n_bases)

    Resource Options:
    --max_cpus      Maximum number of CPUs (int)
                    (default: $params.max_cpus)  
    --max_memory    Maximum memory (memory unit)
                    (default: $params.max_memory)
    --max_time      Maximum time (time unit)
                    (default: $params.max_time)
    
    See here for more info: 
    """.stripIndent()
}

// Show help message
if (params.help) {
  helpMessage()
  exit 0
}

// Define channels from repository files
projectDir = workflow.projectDir

// Define Channels from input
Channel
    .fromPath(params.reference)
    .ifEmpty { exit 1, "Error: Cannot find reference file : ${params.reference}" }
    .set { ch_input_reference }
Channel
    .fromPath(params.annotation)
    .ifEmpty { exit 1, "Error: Cannot find annotation file : ${params.annotation}" }
    .set { ch_input_gtf }

// Define Channels from options
if (!params.read_length && !params.read_length_from && !params.read_length_to){exit 1, "Error: Read length parameters are missing. Either --read_length or --read_length_to and --read_length_from should be provided"}
if (params.read_length_from && !params.read_length_to){exit 1, "Error: Missing --read_length_to when using --read_length_from"}
if (!params.read_length_from && params.read_length_to){exit 1, "Error: Missing --read_length_from when using --read_length_to"}

if (params.read_length) {
  Channel
      .from(params.read_length.split(','))
      .set { ch_read_length }
}

if (params.read_length_from && params.read_length_to) {
  Channel
      .from((params.read_length_from..params.read_length_to).by(params.read_length_by))
      .set { ch_read_length }
}


process STAR {
    tag "$read_length"
    label 'low_memory'
    publishDir "${params.outdir}", mode: 'copy'

    input:
    each file(reference) from ch_input_reference
    each file(gtf) from ch_input_gtf
    val(read_length) from ch_read_length
    
    output:
    file("star_*.tar.gz")

    script:
    ref_name = reference.simpleName

    if (reference.extension == 'gz') {
      ref_file = "${ref_name}.fa"
      unzip_fa_cmd = "gunzip -c ${reference} > ${ref_file}"
    } else {
      ref_file = reference.name
      unzip_fa_cmd = ""
    }

    if (gtf.extension == 'gz') {
      gtf_file = "${gtf.simpleName}.gtf"
      unzip_gtf_cmd = "gunzip -c ${gtf} > ${gtf_file}"
    } else {
      gtf_file = gtf.name
      unzip_gtf_cmd = ""
    }

    """
    ${unzip_fa_cmd}
    ${unzip_gtf_cmd}

    star_version=\$(STAR --version)

    mkdir star_\${star_version}_${ref_name}_$read_length 

    STAR --runMode genomeGenerate \
        --runThreadN $task.cpus \
        --genomeDir star_\${star_version}_${ref_name}_$read_length \
        --genomeFastaFiles ${ref_file} \
        --sjdbGTFfile ${gtf_file} \
        --sjdbOverhang \$(($read_length-1)) \
        --genomeSAindexNbases $params.index_n_bases
    
    tar -czvf star_\${star_version}_${ref_name}_${read_length}.tar.gz star_\${star_version}_${ref_name}_${read_length}

    """
}
