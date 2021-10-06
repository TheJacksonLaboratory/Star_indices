# Star-indices

Workflow for the automation of STAR indices creation from reference data. 

## Requirements

This workflow requires at least XX CPUs and XXGB of memory for optimal execution, however on small input data like test data it can be run with as few as 2 CPUs and 4GB of memory.

## Usage

The typical command for running the pipeline is as follows:

    nextflow run main.nf --reference sample.fa --gtf sample.gtf [Options]

    Inputs Options:
    --reference         Reference file (fasta, gzipped or not)
    --gtf               Annotation file (GTF)

    Read Length Options:
    Use either:
    --read_length       Comma-separated list with read lengths (string)
                        (default: false)
    Or:
    --read_length_from  Start value for read length (int)
                        (default: false)
    --read_length_to    End value for read length (int)
                        (default: false)
    --read_length_by    Value to increment the read length by (int)
                        If not provided, it defaults to an increment of 10.
                        (default: false)

    STAR Options:
    --star_container    Nextflow-compatible container with STAR (string)
                        (default: quay.io/lifebitai/star_indeces:2.7.9a)
    --index_n_bases     Value for STAR genomeSAindexNbases parameter (int)
                        (default: 8)

    Resource Options:
    --max_cpus      Maximum number of CPUs (int)
                    (default: 2)  
    --max_memory    Maximum memory (memory unit)
                    (default: 4 GB)
    --max_time      Maximum time (time unit)
                    (default: 8h)

## Basic run command example

    nextflow run main.nf --reference sample.fa --gtf sample.gtf --read_length_from 50 --read_length_to 150

## Run test locally

    nextflow run main.nf --config conf/test.config

