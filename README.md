# Star-indices

Workflow for the automation of STAR indices creation from reference data. 

## Requirements

This workflow minimal requirements are at least 1 CPUs and 20-60GB of memory (depending on reference genome). For optimal execution providing more resources will speedup the workflow by parallelizing the indexing tasks - 4 cores and 20-60Gb RAM per each read_length you want to run in parallel. For small input data like test data it can be run with as few as 2 CPUs and 4GB of memory.

## Usage

The typical command for running the pipeline is as follows:

    nextflow run main.nf --reference sample.fa --annotation sample.gtf [Options]

    Inputs Options:
    --reference         Reference file (fasta, can be .gz-compressed)
    --annotation        Annotation file (GTF, can be .gz-compressed)

    Read Length Options:
    Use either:
    --read_length       Comma-separated list with read lengths (string)
                        Example: "100,105,110,115,120,150"
                        (default: false)
    Or:
    --read_length_from  Start value for read length (int)
                        (default: false)
    --read_length_to    End value for read length (int)
                        (default: false)
    --read_length_by    Value to increment the read length by (int)
                        If not provided, it defaults to an increment of 10.
                        (default: 10)

    STAR Options:
    --star_container    Nextflow-compatible container with STAR (string)
                        Provided:
                           - quay.io/lifebitai/star_indices:2.7.9a
                           - quay.io/lifebitai/star_indices:2.7.3a
                        (default: quay.io/lifebitai/star_indices:2.7.9a)
    --index_n_bases     Value for STAR genomeSAindexNbases parameter (int)
                        (default: 14)

    Resource Options:
    --max_cpus      Maximum number of CPUs per process (int)
                    (default: 16)  
    --max_memory    Maximum memory per process (memory unit)
                    (default: 100.GB)
    --max_time      Maximum time per process (time unit)
                    (default: 72.h)

## Basic run command example

    nextflow run main.nf --reference sample.fa --annotation sample.gtf --read_length_from 50 --read_length_to 150

## Run test locally

    nextflow run main.nf --config conf/test.config

