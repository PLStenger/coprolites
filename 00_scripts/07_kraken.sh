#!/bin/bash
#SBATCH --job-name=07_kraken
#SBATCH --ntasks=1
#SBATCH -p smp
#SBATCH --mem=1000G
#SBATCH --mail-user=pierrelouis.stenger@gmail.com
#SBATCH --mail-type=ALL 
#SBATCH --error="/home/plstenge/coprolites/00_scripts/07_kraken.err"
#SBATCH --output="/home/plstenge/coprolites/00_scripts/07_kraken.out"

module load conda/4.12.0
source ~/.bashrc
conda activate kraken2


### Download Kraken database from https://benlangmead.github.io/aws-indexes/k2
### k2_core_nt_20250609.tar.gz

# Répertoires
FASTP_DIR="/home/plstenge/coprolites/06_fastp"
KRAKEN2_DB="/home/plstenge/k2_core_nt_20250609"
OUT_DIR="/home/plstenge/coprolites/07_kraken2"
THREADS=36

mkdir -p "$OUT_DIR"

# 1. Analyse des merged (single-end)
echo "=== Analyse Kraken2 des merged (single-end) ==="
for MERGED in "$FASTP_DIR"/*_fastp_merged.fastq.gz
do
    SAMPLE=$(basename "$MERGED" _fastp_merged.fastq.gz)
    OUT_KRAKEN="$OUT_DIR/${SAMPLE}_merged.kraken"
    OUT_REPORT="$OUT_DIR/${SAMPLE}_merged.report"

    kraken2 --conf 0.2 --db "$KRAKEN2_DB" --threads $THREADS \
        --output "$OUT_KRAKEN" --report "$OUT_REPORT" "$MERGED"

    echo "Terminé : $SAMPLE (merged)"
done

# 2. Analyse des unmerged (paired-end)
echo "=== Analyse Kraken2 des unmerged (paired-end) ==="
for R1 in "$FASTP_DIR"/*_fastp_R1.fastq.gz
do
    SAMPLE=$(basename "$R1" _fastp_R1.fastq.gz)
    R2="$FASTP_DIR/${SAMPLE}_fastp_R2.fastq.gz"

    # Ne lance que si R2 existe
    if [[ -f "$R2" ]]; then
        OUT_KRAKEN="$OUT_DIR/${SAMPLE}_unmerged.kraken"
        OUT_REPORT="$OUT_DIR/${SAMPLE}_unmerged.report"

        kraken2 --conf 0.2 --paired --db "$KRAKEN2_DB" --threads $THREADS \
            --output "$OUT_KRAKEN" --report "$OUT_REPORT" "$R1" "$R2"

        echo "Terminé : $SAMPLE (unmerged)"
    fi
done

echo "Analyse Kraken2 terminée pour tous les échantillons."
