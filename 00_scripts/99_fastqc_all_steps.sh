#!/bin/bash
#SBATCH --job-name=quality_check_by_FastQC
##SBATCH --time=24:00
#SBATCH --ntasks=1
#SBATCH -p smp
##SBATCH --nodelist=gdecnode02
#SBATCH --mem=250G
##SBATCH -c 32
#SBATCH --mail-user=pierrelouis.stenger@gmail.com
#SBATCH --mail-type=ALL 
#SBATCH --error="/home/plstenge/coprolites/00_scripts/99_fastqc_all_steps.err"
#SBATCH --output="/home/plstenge/coprolites/00_scripts/99_fastqc_all_steps.out"

# installing FastQC from https://www.bioinformatics.babraham.ac.uk/projects/download.html
# FastQC v0.11.9 (Mac DMG image)

# Correct tool citation : Andrews, S. (2010). FastQC: a quality control tool for high throughput sequence data.

###############################################################################################################################
######## Step 1: bbduk #######
###############################################################################################################################
#
#WORKING_DIRECTORY=/home/plstenge/coprolites/03_bbduk
#OUTPUT=/home/plstenge/coprolites/99_02_quality_check_bbduk
#
## Make the directory (mkdir) only if not existe already(-p)
#mkdir -p $OUTPUT
#
#module load conda/4.12.0
#source ~/.bashrc
#conda activate fastqc
#
#cd $WORKING_DIRECTORY
#
#for FILE in $(ls $WORKING_DIRECTORY/*.fastq.gz)
#do
#      fastqc $FILE -o $OUTPUT
#done ;
#
#conda deactivate
#
#module load conda/4.12.0
#source ~/.bashrc
#conda activate multiqc
#
## Run multiqc for quality summary
#
#multiqc $OUTPUT
#
#conda deactivate
#
#
###############################################################################################################################
######## Step 2: fastuniq #######
###############################################################################################################################
#
#WORKING_DIRECTORY=/home/plstenge/coprolites/04_fastuniq
#OUTPUT=/home/plstenge/coprolites/99_03_quality_check_fastuniq
#
## Make the directory (mkdir) only if not existe already(-p)
#mkdir -p $OUTPUT
#
#module load conda/4.12.0
#source ~/.bashrc
#conda activate fastqc
#
#cd $WORKING_DIRECTORY
#
#for FILE in $(ls $WORKING_DIRECTORY/*.fastq.gz)
#do
#      fastqc $FILE -o $OUTPUT
#done ;
#
#conda deactivate
#
#module load conda/4.12.0
#source ~/.bashrc
#conda activate multiqc
#
## Run multiqc for quality summary
#
#multiqc $OUTPUT
#
#conda deactivate
#
#
###############################################################################################################################
######## Step 3: clumpify #######
###############################################################################################################################
#
#WORKING_DIRECTORY=/home/plstenge/coprolites/05_clumpify
#OUTPUT=/home/plstenge/coprolites/99_04_quality_check_clumpify
#
## Make the directory (mkdir) only if not existe already(-p)
#mkdir -p $OUTPUT
#
#module load conda/4.12.0
#source ~/.bashrc
#conda activate fastqc
#
#cd $WORKING_DIRECTORY
#
#for FILE in $(ls $WORKING_DIRECTORY/*.fastq.gz)
#do
#      fastqc $FILE -o $OUTPUT
#done ;
#
#conda deactivate
#
#module load conda/4.12.0
#source ~/.bashrc
#conda activate multiqc
#
## Run multiqc for quality summary
#
#multiqc $OUTPUT
#
#conda deactivate

##############################################################################################################################
####### Step 4: fastp #######
##############################################################################################################################

WORKING_DIRECTORY=/home/plstenge/coprolites/06_fastp
OUTPUT=/home/plstenge/coprolites/99_05_quality_check_fastp

# Make the directory (mkdir) only if not existe already(-p)
mkdir -p $OUTPUT

module load conda/4.12.0
source ~/.bashrc
conda activate fastqc

cd $WORKING_DIRECTORY

for FILE in $(ls $WORKING_DIRECTORY/*.fastq.gz)
do
      fastqc $FILE -o $OUTPUT
done ;

conda deactivate

module load conda/4.12.0
source ~/.bashrc
conda activate multiqc

# Run multiqc for quality summary

multiqc $OUTPUT

conda deactivate
