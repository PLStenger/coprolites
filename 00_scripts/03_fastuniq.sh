#!/bin/bash
#SBATCH --job-name=03_fastuniq
#SBATCH --ntasks=1
#SBATCH -p smp
#SBATCH --mem=250G
#SBATCH --mail-user=pierrelouis.stenger@gmail.com
#SBATCH --mail-type=ALL 
#SBATCH --error="/home/plstenge/coprolites/00_scripts/03_fastuniq.err"
#SBATCH --output="/home/plstenge/coprolites/00_scripts/03_fastuniq.out"

INPUT="/home/plstenge/coprolites/03_bbduk"
OUTPUT="/home/plstenge/coprolites/04_fastuniq"

mkdir -p "$OUTPUT"

module load conda/4.12.0
source ~/.bashrc
conda activate fastuniq

cd "$INPUT" || exit 1

TMP="/tmp/fastuniq_tmp"
mkdir -p "$TMP"

for R1_gz in clean_*_R1.fastq.gz; do
    base=$(echo "$R1_gz" | sed 's/_R1\.fastq\.gz//')
    R2_gz="${base}_R2.fastq.gz"

    if [[ -f "$R2_gz" ]]; then
        echo "Traitement de la paire: $base"

        R1_tmp="${TMP}/${base}_R1.fastq"
        R2_tmp="${TMP}/${base}_R2.fastq"
        listfile="${TMP}/${base}.list"

        zcat "$INPUT/$R1_gz" > "$R1_tmp"
        zcat "$INPUT/$R2_gz" > "$R2_tmp"

        echo -e "$R1_tmp\n$R2_tmp" > "$listfile"

        fastuniq -i "$listfile" -t q \
            -o "${OUTPUT}/${base}_dedup_R1.fastq" \
            -p "${OUTPUT}/${base}_dedup_R2.fastq"

        rm -f "$R1_tmp" "$R2_tmp" "$listfile"

    else
        echo "ATTENTION: fichier R2 manquant pour $base"
    fi
done

echo "Terminé."
