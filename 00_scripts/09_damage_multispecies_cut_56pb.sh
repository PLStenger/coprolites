#!/bin/bash
#SBATCH --job-name=09_damage_multispecies_cut_56pb
#SBATCH --ntasks=1
#SBATCH -p smp
#SBATCH --mem=1000G
#SBATCH --mail-user=pierrelouis.stenger@gmail.com
#SBATCH --mail-type=ALL 
#SBATCH --error="/home/plstenge/coprolites/00_scripts/09_damage_multispecies_cut_56pb.err"
#SBATCH --output="/home/plstenge/coprolites/00_scripts/09_damage_multispecies_cut_56pb.out"

module load conda/4.12.0
source ~/.bashrc
conda activate mapdamage_py39

KRAKEN_DIR="/home/plstenge/coprolites/07_kraken2"
FASTQ_DIR="/home/plstenge/coprolites/06_fastp"
DAMAGE_BASE="/home/plstenge/coprolites/08_damage"
LOGFILE="/home/plstenge/coprolites/00_scripts/09_HOPS_homemade_$(date +%Y%m%d_%H%M%S).txt"
MAPPING_SUMMARY="/home/plstenge/coprolites/08_damage/mapping_rates.txt"

echo "Script started at $(date)" | tee -a "$LOGFILE"
bwa index /home/plstenge/genomes/Ovis_aries.ARS-UI_Ramb_v3.0.dna.toplevel.fa
#bwa index /home/plstenge/genomes/Capra_hircus.ARS1.dna.toplevel.fa
#bwa index /home/plstenge/genomes/Alnus_glutinosa_genome_assembly_dhAlnGlut1.fa
#bwa index /home/plstenge/genomes/Corylus_avellana_CavTom2PMs_1_0.fasta

declare -A TAXONS=(
    ["Ovis_aries"]="9940:/home/plstenge/genomes/Ovis_aries.ARS-UI_Ramb_v3.0.dna.toplevel.fa"
#    ["Capra_hircus"]="9925:/home/plstenge/genomes/Capra_hircus.ARS1.dna.toplevel.fa"
#    ["Alnus_glutinosa"]="3517:/home/plstenge/genomes/Alnus_glutinosa_genome_assembly_dhAlnGlut1.fa"
#    ["Corylus_avellana"]="13451:/home/plstenge/genomes/Corylus_avellana_CavTom2PMs_1_0.fasta"
)
shopt -s nullglob

echo -e "Sample\tGroup\tUnmerged_MappingRate\tMerged_MappingRate" > "$MAPPING_SUMMARY"

for KRAKEN_FILE in "$KRAKEN_DIR"/*.kraken; do
    KRAKEN_BASE=$(basename "$KRAKEN_FILE" .kraken)
    echo -e "\n==== Processing file: $KRAKEN_FILE ====" | tee -a "$LOGFILE"
    PREFIX=$(echo "$KRAKEN_BASE" | sed -E 's/_dedup_clumpify_(un)?merged$//')
    echo "Prefix for FASTQ search: $PREFIX" | tee -a "$LOGFILE"
    R1_FILES=("${FASTQ_DIR}/${PREFIX}_dedup_clumpify_fastp_R1.fastq"*)
    R2_FILES=("${FASTQ_DIR}/${PREFIX}_dedup_clumpify_fastp_R2.fastq"*)
    MERGED_FILES=("${FASTQ_DIR}/${PREFIX}_dedup_clumpify_fastp_merged.fastq"*)
    R1_FILE="${R1_FILES[0]:-}"
    R2_FILE="${R2_FILES[0]:-}"
    MERGED_FILE="${MERGED_FILES[0]:-}"

    for GROUP in "${!TAXONS[@]}"; do
        TAX_ID="${TAXONS[$GROUP]%:*}"
        REF_FASTA="${TAXONS[$GROUP]#*:}"
        DAMAGE_DIR="${DAMAGE_BASE}/${GROUP}"
        mkdir -p "$DAMAGE_DIR"
        
        OUT_R1="${DAMAGE_DIR}/${KRAKEN_BASE}_${GROUP}_R1.fastq"
        OUT_R2="${DAMAGE_DIR}/${KRAKEN_BASE}_${GROUP}_R2.fastq"
        OUT_MERGED="${DAMAGE_DIR}/${KRAKEN_BASE}_${GROUP}_merged.fastq"

        # Extraction Kraken
        if [[ -n "$R1_FILE" && -n "$R2_FILE" ]]; then
            python3 /home/plstenge/KrakenTools/extract_kraken_reads.py \
                -k "$KRAKEN_FILE" -s "$R1_FILE" -s2 "$R2_FILE" -t "$TAX_ID" \
                -o "$OUT_R1" -o2 "$OUT_R2" --fastq-output 2>>"$LOGFILE"

            # Trimming si "Aviti"
            if [[ "$KRAKEN_BASE" =~ Aviti ]]; then
                seqtk trimfq -b 56 "$OUT_R1" > "${OUT_R1}.trim"
                seqtk trimfq -b 56 "$OUT_R2" > "${OUT_R2}.trim"
                OUT_R1="${OUT_R1}.trim"
                OUT_R2="${OUT_R2}.trim"
            fi

            if [[ -f "$OUT_R1" && -f "$OUT_R2" ]]; then
                bwa aln -n 0.08 -l 24 -k 2 -q 20 -t 4 "$REF_FASTA" "$OUT_R1" > "${DAMAGE_DIR}/${KRAKEN_BASE}_${GROUP}_R1.sai" 2>>"$LOGFILE"
                bwa aln -n 0.08 -l 24 -k 2 -q 20 -t 4 "$REF_FASTA" "$OUT_R2" > "${DAMAGE_DIR}/${KRAKEN_BASE}_${GROUP}_R2.sai" 2>>"$LOGFILE"
                bwa sampe "$REF_FASTA" "${DAMAGE_DIR}/${KRAKEN_BASE}_${GROUP}_R1.sai" "${DAMAGE_DIR}/${KRAKEN_BASE}_${GROUP}_R2.sai" "$OUT_R1" "$OUT_R2" > "${DAMAGE_DIR}/${KRAKEN_BASE}_${GROUP}.sam" 2>>"$LOGFILE"
                samtools view -bS "${DAMAGE_DIR}/${KRAKEN_BASE}_${GROUP}.sam" > "${DAMAGE_DIR}/${KRAKEN_BASE}_${GROUP}.bam" 2>>"$LOGFILE"
                samtools sort -o "${DAMAGE_DIR}/${KRAKEN_BASE}_${GROUP}.sorted.bam" "${DAMAGE_DIR}/${KRAKEN_BASE}_${GROUP}.bam" 2>>"$LOGFILE"
                samtools index "${DAMAGE_DIR}/${KRAKEN_BASE}_${GROUP}.sorted.bam" 2>>"$LOGFILE"
                rm -f ${DAMAGE_DIR}/${KRAKEN_BASE}_${GROUP}_R1.sai ${DAMAGE_DIR}/${KRAKEN_BASE}_${GROUP}_R2.sai ${DAMAGE_DIR}/${KRAKEN_BASE}_${GROUP}.sam ${DAMAGE_DIR}/${KRAKEN_BASE}_${GROUP}.bam 2>>"$LOGFILE"
                mapDamage -i "${DAMAGE_DIR}/${KRAKEN_BASE}_${GROUP}.sorted.bam" -r "$REF_FASTA" --folder "${DAMAGE_DIR}/${KRAKEN_BASE}_${GROUP}_mapDamage_unmerged" 2>>"$LOGFILE"
                # Taux de mapping pour les unmerged
                UNMERGED_MAPRATE=$(samtools flagstat "${DAMAGE_DIR}/${KRAKEN_BASE}_${GROUP}.sorted.bam" | grep "mapped (" | head -1 | awk '{print $5}')
            else
                UNMERGED_MAPRATE="NA"
            fi
        else
            UNMERGED_MAPRATE="NA"
        fi

        # Merged
        if [[ -n "$MERGED_FILE" ]]; then
            python3 /home/plstenge/KrakenTools/extract_kraken_reads.py \
                -k "$KRAKEN_FILE" -s "$MERGED_FILE" -t "$TAX_ID" -o "$OUT_MERGED" --fastq-output 2>>"$LOGFILE"

            # Trimming si "Aviti"
            if [[ "$KRAKEN_BASE" =~ Aviti ]]; then
                seqtk trimfq -b 56 "$OUT_MERGED" > "${OUT_MERGED}.trim"
                OUT_MERGED="${OUT_MERGED}.trim"
            fi

            if [[ -f "$OUT_MERGED" ]]; then
                bwa aln -n 0.08 -l 24 -k 2 -q 20 -t 4 "$REF_FASTA" "$OUT_MERGED" > "${DAMAGE_DIR}/${KRAKEN_BASE}_${GROUP}_merged.sai" 2>>"$LOGFILE"
                bwa samse "$REF_FASTA" "${DAMAGE_DIR}/${KRAKEN_BASE}_${GROUP}_merged.sai" "$OUT_MERGED" > "${DAMAGE_DIR}/${KRAKEN_BASE}_${GROUP}_merged.sam" 2>>"$LOGFILE"
                samtools view -bS "${DAMAGE_DIR}/${KRAKEN_BASE}_${GROUP}_merged.sam" > "${DAMAGE_DIR}/${KRAKEN_BASE}_${GROUP}_merged.bam" 2>>"$LOGFILE"
                samtools sort -o "${DAMAGE_DIR}/${KRAKEN_BASE}_${GROUP}_merged.sorted.bam" "${DAMAGE_DIR}/${KRAKEN_BASE}_${GROUP}_merged.bam" 2>>"$LOGFILE"
                samtools index "${DAMAGE_DIR}/${KRAKEN_BASE}_${GROUP}_merged.sorted.bam" 2>>"$LOGFILE"
                rm -f ${DAMAGE_DIR}/${KRAKEN_BASE}_${GROUP}_merged.sai ${DAMAGE_DIR}/${KRAKEN_BASE}_${GROUP}_merged.sam ${DAMAGE_DIR}/${KRAKEN_BASE}_${GROUP}_merged.bam 2>>"$LOGFILE"
                mapDamage -i "${DAMAGE_DIR}/${KRAKEN_BASE}_${GROUP}_merged.sorted.bam" -r "$REF_FASTA" --folder "${DAMAGE_DIR}/${KRAKEN_BASE}_${GROUP}_mapDamage_merged" 2>>"$LOGFILE"
                # Taux de mapping pour les merged
                MERGED_MAPRATE=$(samtools flagstat "${DAMAGE_DIR}/${KRAKEN_BASE}_${GROUP}_merged.sorted.bam" | grep "mapped (" | head -1 | awk '{print $5}')
            else
                MERGED_MAPRATE="NA"
            fi
        else
            MERGED_MAPRATE="NA"
        fi

        # Résumé dans un fichier table
        echo -e "${KRAKEN_BASE}\t${GROUP}\t${UNMERGED_MAPRATE}\t${MERGED_MAPRATE}" >> "$MAPPING_SUMMARY"
    done
done

echo "Finished at $(date)" | tee -a "$LOGFILE"
