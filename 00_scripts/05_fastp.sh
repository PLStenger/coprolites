#!/bin/bash
#SBATCH --job-name=05_fastp
#SBATCH --ntasks=1
#SBATCH -p smp
#SBATCH --mem=250G
#SBATCH --mail-user=pierrelouis.stenger@gmail.com
#SBATCH --mail-type=ALL 
#SBATCH --error="/home/plstenge/coprolites/00_scripts/05_fastp.err"
#SBATCH --output="/home/plstenge/coprolites/00_scripts/05_fastp.out"

# Dossier contenant les fichiers fastq.gz
INPUT_DIR="/home/plstenge/coprolites/05_clumpify"
OUTPUT_DIR="/home/plstenge/coprolites/06_fastp"
LOG_DIR="/home/plstenge/coprolites/00_scripts/05_fastp_out/"

# Création des dossiers de sortie si besoin
mkdir -p "$OUTPUT_DIR"
mkdir -p "$LOG_DIR"

module load conda/4.12.0
source ~/.bashrc
conda activate fastp

# Boucle sur tous les fichiers R1
for R1 in "$INPUT_DIR"/*_R1.fastq.gz
do
    # Trouver le fichier R2 correspondant
    R2="${R1/_R1.fastq.gz/_R2.fastq.gz}"
    
    # Nom de base pour les sorties
    BASENAME=$(basename "$R1" _R1.fastq.gz)
    
    # Fichiers de sortie
    OUT_R1="$OUTPUT_DIR/${BASENAME}_fastp_R1.fastq.gz"
    OUT_R2="$OUTPUT_DIR/${BASENAME}_fastp_R2.fastq.gz"
    MERGED="$OUTPUT_DIR/${BASENAME}_fastp_merged.fastq.gz"
    HTML="$LOG_DIR/${BASENAME}_fastp.html"
    JSON="$LOG_DIR/${BASENAME}_fastp.json"
    
    # Lancer fastp
    fastp \
        --in1 "$R1" --in2 "$R2" \
        --out1 "$OUT_R1" --out2 "$OUT_R2" \
        --merged_out "$MERGED" \
        --length_required 30 \
        --cut_front --cut_tail \
        --cut_window_size 4 \
        --cut_mean_quality 10 \
        --n_base_limit 5 \
        --unqualified_percent_limit 40 \
        --complexity_threshold 30 \
        --qualified_quality_phred 15 \
        --low_complexity_filter \
        --trim_poly_x \
        --poly_x_min_len 10 \
        --merge --correction \
        --overlap_len_require 10 \
        --overlap_diff_limit 5 \
        --overlap_diff_percent_limit 20 \
        --html "$HTML" \
        --json "$JSON" \
        --adapter_sequence AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
        --adapter_sequence_r2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
        --detect_adapter_for_pe \
        --thread 4

    echo "Traitement terminé pour $BASENAME"
done

echo "Tous les fichiers ont été traités."
