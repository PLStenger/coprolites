#!/bin/bash
#SBATCH --job-name=09_HOPS_reduced_02
#SBATCH --ntasks=1
#SBATCH -p smp
#SBATCH --mem=1000G
#SBATCH --mail-user=pierrelouis.stenger@gmail.com
#SBATCH --mail-type=ALL 
#SBATCH --error="/home/plstenge/coprolites/00_scripts/09_HOPS_reduced_02.err"
#SBATCH --output="/home/plstenge/coprolites/00_scripts/09_HOPS_reduced_02.out"

set -euo pipefail

KR_DIR="/home/plstenge/coprolites/07_kraken2"
OUT_DIR="/home/plstenge/custom_nt_hopped"

module load conda/4.12.0
source ~/.bashrc
conda activate malt

mkdir -p $OUT_DIR
cd $OUT_DIR

# Pour augmenter la memoire de l'outil (bridée à 64Go):
# cd /home/plstenge/miniconda3/envs/malt
# find . -name "*.vmoptions"
# Modifier le XmX

export JAVA_TOOL_OPTIONS="-Xmx1000G"

# -----------------------------------------------------------
#     --- Construction de l’index MALT personnalisé ---
# -----------------------------------------------------------
malt-build -a2t final_accession2taxid.txt \
        --index HOPS_custom_nt_filtered \
        --input custom_nt.fasta \
        --threads 16 \
        --memoryMode load
