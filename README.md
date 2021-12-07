# Star-indices

Workflow for the automation of STAR indices creation from reference data. 

## Hardware requirements

This workflow minimal requirements highly depend on the reference organism genome size, especially for RAM usage. For CPUs at least 1 CPU is enough to run the pipeline, however it does highly benefit from parallelisation on up to 16-32 CPUs. Below are the minmal and optimal resource requirements needed to run the pipeline locally and on the cloud.

### Minimal requirement to run locally

| Organism |Genome fasta size, Gb | Max average CPU utilisation observed | Minimal RAM | Estimated runtime |
|---|---|---|---|---|
| Yeast      | 0.01 | 1 | 4 Gb  | 0:04:00 |
| Drosophila | 0.15 | 2 | 7 Gb  | 0:10:00 |
| Mouse      | 2.78 | 4 | 29 Gb | 1:20:00 |
| Human      | 3.15 | 8 | 34 Gb | 1:30:00 |

Note that for higher overhang sizes and higher number of CPUs involved RAM usage can go up to 10-20% higher than stated in the table above.

### Optimal requirements to run on the Google Cloud
On Google Cloudit it is not possible to request number of CPUs and gigabytes of RAM very different from proportion 1 CPU - 4 Gb RAM (such as 1 CPU - 32 GB). Therefore optimal resource usage on cloud must include more CPUs for higher RAM usage even if pipeline cannot take advantage of all CPUs.

| Organism |Genome fasta size, Gb | Optimal CPUs | Optimal RAM | Optimal runtime |
|---|---|---|---|---|
| Yeast      | 0.01 | 2 | 4 Gb | 0:04:00 |
| Drosophila | 0.15 | 4 | 8 Gb | 0:10:00 |
| Mouse      | 2.78 | 8 | 32 Gb | 1:20:00 |
| Human      | 3.15 |16 | 64 Gb | 1:30:00 |

Supporting information: benchmarking [google sheet](https://docs.google.com/spreadsheets/d/1s6D8XDZQFxa-ac8IcGjHe_iJsdq8RwZqzAWxnjTNkfc/edit?usp=sharing).

For small input data like test data it can be run with as few as 2 CPUs and 4GB of memory.

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

