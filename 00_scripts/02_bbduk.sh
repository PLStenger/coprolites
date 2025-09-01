#!/bin/bash
#SBATCH --job-name=02_bbduk
#SBATCH --ntasks=1
#SBATCH -p smp
#SBATCH --mem=250G
#SBATCH --mail-user=pierrelouis.stenger@gmail.com
#SBATCH --mail-type=ALL 
#SBATCH --error="/home/plstenge/coprolites/00_scripts/02_bbduk.err"
#SBATCH --output="/home/plstenge/coprolites/00_scripts/02_bbduk.out"

WORKING_DIRECTORY="/home/plstenge/coprolites/01_raw_data"
OUTPUT="/home/plstenge/coprolites/03_bbduk"

mkdir -p "$OUTPUT"

module load conda/4.12.0
source ~/.bashrc
conda activate bbduk

PHIX=/home/plstenge/bbmap/resources/phix174_ill.ref.fa.gz
BBDUK=/home/plstenge/bbmap/bbduk.sh

cd "$WORKING_DIRECTORY"

# Boucle 1: Démultiplexage et trimming/collapse pour chaque échantillon
for r1_file in *_R1.fastq.gz; do
    r2_file="${r1_file/_R1/_R2}"
    [[ ! -f "$r2_file" ]] && { echo "ERREUR: Fichier R2 manquant pour $r1_file" >&2; continue; }
    base_name="${r1_file%%_R1.fastq.gz}"
  
   $BBDUK -Xmx4g \
    in1="$r1_file" \
    in2="$r2_file" \
    out1=$OUTPUT/"clean_${r1_file}" \
    out2=$OUTPUT/"clean_${r2_file}" \
    ref=$PHIX \
    ktrim=rl \
    k=23 \
    mink=11 \
    hdist=1 \
    tpe \
    tbo \
    minlen=25 \
    qtrim=r \
    trimq=20 \
    stats="${r1_file%.*}_stats.txt"
    
done
