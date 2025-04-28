#!/bin/bash
#SBATCH -p short
#SBATCH --job-name=mouse_wt_long_dox_timecourse
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=reow9695@colorado.edu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=6gb
#SBATCH --time=10:00:00
#SBATCH --output=nextflow.out
#SBATCH --error=netflow.err

pwd; hostname; date
echo "yipppeeeeee lesson 4 is a go"

module load singularity/3.1.1
module load openjdk/21.0.1

nextflow run nf-core/rnaseq -r 3.14.0 \
-resume \
-profile singularity \
--input samplesheet.csv \
--outdir /scratch/Shares/rinnclass/MASTER_CLASS/STUDENTS/reanna/MASTER_CLASS/lessons/04_RNAseq_Dox/01_Mouse_dox_wt/dox_long/dox_long_output \
--reads /scratch/Shares/rinnclass/MASTER_CLASS/DATA/mouse_wt_long_timecourse/*{_R1,_R2}.fastq.gz \
--fasta /scratch/Shares/rinnclass/MASTER_CLASS/GENOMES/M25/GRCm38.p6.genome.fa \
--gtf /scratch/Shares/rinnclass/MASTER_CLASS/GENOMES/M25/gencode.vM25.annotation.gtf \
--psuedo_aligner salmon \
--gencode \
--email reow9695@colorado.edu
-c nextflox.config
