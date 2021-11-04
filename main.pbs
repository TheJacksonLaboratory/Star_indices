#!/bin/bash
#SBATCH -o splicing.%j.out
#SBATCH -e splicing.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=$USER@jax.org
#SBATCH --mem=35000
#SBATCH --cpus-per-task=4
#SBATCH -p compute
#SBATCH -q batch
#SBATCH -t 1-00:00:00 

cd $SLURM_SUBMIT_DIR
date;hostname;pwd

module load singularity

curl -fsSL get.nextflow.io | bash

./nextflow run /projects/anczukow-lab/star_index_pipeline/Star_indices/main.nf \
	--outdir ${SLURM_SUBMIT_DIR} \
	-config NF_Star_Index.config \
	-profile sumner -resume \
        -with-report report.html \
        -with-timeline timeline.html 
