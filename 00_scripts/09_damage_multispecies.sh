#!/bin/bash
#SBATCH --job-name=09_damage_multispecies
#SBATCH --ntasks=1
#SBATCH -p smp
#SBATCH --mem=1000G
#SBATCH --mail-user=pierrelouis.stenger@gmail.com
#SBATCH --mail-type=ALL 
#SBATCH --error="/home/plstenge/coprolites/00_scripts/09_damage_multispecies.err"
#SBATCH --output="/home/plstenge/coprolites/00_scripts/09_damage_multispecies.out"

module load conda/4.12.0
source ~/.bashrc
conda activate mapdamage_py39

KRAKEN_DIR="/home/plstenge/coprolites/07_kraken2"
FASTQ_DIR="/home/plstenge/coprolites/06_fastp"
DAMAGE_BASE="/home/plstenge/coprolites/08_damage"
LOGFILE="/home/plstenge/coprolites/00_scripts/09_HOPS_homemade_$(date +%Y%m%d_%H%M%S).txt"
MAPPING_INFO="/home/plstenge/coprolites/00_scripts/mapping_bwa_info.txt"

echo "Script started at $(date)" | tee -a "$LOGFILE"

# Initialiser le fichier de mapping info avec les en-têtes
echo -e "Sample\tSpecies\tType\tTotal_Reads\tMapped_Reads\tMapping_Rate(%)" > "$MAPPING_INFO"

#bwa index /home/plstenge/genomes/Ovis_aries.ARS-UI_Ramb_v3.0.dna.toplevel.fa
#bwa index /home/plstenge/genomes/Capra_hircus.ARS1.dna.toplevel.fa
#bwa index /home/plstenge/genomes/Alnus_glutinosa_genome_assembly_dhAlnGlut1.fa
#bwa index /home/plstenge/genomes/Corylus_avellana_CavTom2PMs_1_0.fasta
bwa index /home/plstenge/genomes/Phragmites_australis_GCA_040373225.1.fasta

declare -A TAXONS=(
 #   ["Ovis_aries"]="9940:/home/plstenge/genomes/Ovis_aries.ARS-UI_Ramb_v3.0.dna.toplevel.fa"
 #   ["Capra_hircus"]="9925:/home/plstenge/genomes/Capra_hircus.ARS1.dna.toplevel.fa"
  #  ["Alnus_glutinosa"]="3517:/home/plstenge/genomes/Alnus_glutinosa_genome_assembly_dhAlnGlut1.fa"
    ["Phragmites_australis"]="29695:/home/plstenge/genomes/Phragmites_australis_GCA_040373225.1.fasta"
  #  ["Corylus_avellana"]="13451:/home/plstenge/genomes/Corylus_avellana_CavTom2PMs_1_0.fasta"
)

# Fonction pour calculer le taux de mapping
calculate_mapping_rate() {
    local bam_file="$1"
    local sample_name="$2"
    local species="$3"
    local type="$4"
    
    if [[ -f "$bam_file" ]]; then
        # Compter le nombre total de reads
        local total_reads=$(samtools view -c "$bam_file")
        
        # Compter le nombre de reads mappés (flag != 4)
        local mapped_reads=$(samtools view -c -F 4 "$bam_file")
        
        # Calculer le taux de mapping
        local mapping_rate=0
        if [[ $total_reads -gt 0 ]]; then
            mapping_rate=$(echo "scale=2; $mapped_reads * 100 / $total_reads" | bc)
        fi
        
        # Écrire dans le fichier de résultats
        echo -e "${sample_name}\t${species}\t${type}\t${total_reads}\t${mapped_reads}\t${mapping_rate}" >> "$MAPPING_INFO"
        
        echo "Mapping stats for ${sample_name}_${species}_${type}: ${mapped_reads}/${total_reads} (${mapping_rate}%)" | tee -a "$LOGFILE"
    fi
}

shopt -s nullglob
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
        
        # Unmerged pairs
        if [[ -n "$R1_FILE" && -n "$R2_FILE" ]]; then
            python3 /home/plstenge/KrakenTools/extract_kraken_reads.py \
                -k "$KRAKEN_FILE" -s "$R1_FILE" -s2 "$R2_FILE" -t "$TAX_ID" \
                -o "$OUT_R1" -o2 "$OUT_R2" --fastq-output 2>>"$LOGFILE"
            if [[ -f "$OUT_R1" && -f "$OUT_R2" ]]; then
                bwa aln -n 0.08 -l 24 -k 2 -q 20 -t 4 "$REF_FASTA" "$OUT_R1" > "${DAMAGE_DIR}/${KRAKEN_BASE}_${GROUP}_R1.sai" 2>>"$LOGFILE"
                bwa aln -n 0.08 -l 24 -k 2 -q 20 -t 4 "$REF_FASTA" "$OUT_R2" > "${DAMAGE_DIR}/${KRAKEN_BASE}_${GROUP}_R2.sai" 2>>"$LOGFILE"
                bwa sampe "$REF_FASTA" "${DAMAGE_DIR}/${KRAKEN_BASE}_${GROUP}_R1.sai" "${DAMAGE_DIR}/${KRAKEN_BASE}_${GROUP}_R2.sai" "$OUT_R1" "$OUT_R2" > "${DAMAGE_DIR}/${KRAKEN_BASE}_${GROUP}.sam" 2>>"$LOGFILE"
                samtools view -bS "${DAMAGE_DIR}/${KRAKEN_BASE}_${GROUP}.sam" > "${DAMAGE_DIR}/${KRAKEN_BASE}_${GROUP}.bam" 2>>"$LOGFILE"
                samtools sort -o "${DAMAGE_DIR}/${KRAKEN_BASE}_${GROUP}.sorted.bam" "${DAMAGE_DIR}/${KRAKEN_BASE}_${GROUP}.bam" 2>>"$LOGFILE"
                samtools index "${DAMAGE_DIR}/${KRAKEN_BASE}_${GROUP}.sorted.bam" 2>>"$LOGFILE"
                
                # Calculer le taux de mapping pour les reads unmerged
                calculate_mapping_rate "${DAMAGE_DIR}/${KRAKEN_BASE}_${GROUP}.sorted.bam" "$KRAKEN_BASE" "$GROUP" "unmerged"
                
                rm -f ${DAMAGE_DIR}/${KRAKEN_BASE}_${GROUP}_R1.sai ${DAMAGE_DIR}/${KRAKEN_BASE}_${GROUP}_R2.sai ${DAMAGE_DIR}/${KRAKEN_BASE}_${GROUP}.sam ${DAMAGE_DIR}/${KRAKEN_BASE}_${GROUP}.bam 2>>"$LOGFILE"
                mapDamage -i "${DAMAGE_DIR}/${KRAKEN_BASE}_${GROUP}.sorted.bam" -r "$REF_FASTA" --folder "${DAMAGE_DIR}/${KRAKEN_BASE}_${GROUP}_mapDamage_unmerged" 2>>"$LOGFILE"
            fi
        fi
        # Merged
        if [[ -n "$MERGED_FILE" ]]; then
            python3 /home/plstenge/KrakenTools/extract_kraken_reads.py \
                -k "$KRAKEN_FILE" -s "$MERGED_FILE" -t "$TAX_ID" -o "$OUT_MERGED" --fastq-output 2>>"$LOGFILE"
            if [[ -f "$OUT_MERGED" ]]; then
                bwa aln -n 0.08 -l 24 -k 2 -q 20 -t 4 "$REF_FASTA" "$OUT_MERGED" > "${DAMAGE_DIR}/${KRAKEN_BASE}_${GROUP}_merged.sai" 2>>"$LOGFILE"
                bwa samse "$REF_FASTA" "${DAMAGE_DIR}/${KRAKEN_BASE}_${GROUP}_merged.sai" "$OUT_MERGED" > "${DAMAGE_DIR}/${KRAKEN_BASE}_${GROUP}_merged.sam" 2>>"$LOGFILE"
                samtools view -bS "${DAMAGE_DIR}/${KRAKEN_BASE}_${GROUP}_merged.sam" > "${DAMAGE_DIR}/${KRAKEN_BASE}_${GROUP}_merged.bam" 2>>"$LOGFILE"
                samtools sort -o "${DAMAGE_DIR}/${KRAKEN_BASE}_${GROUP}_merged.sorted.bam" "${DAMAGE_DIR}/${KRAKEN_BASE}_${GROUP}_merged.bam" 2>>"$LOGFILE"
                samtools index "${DAMAGE_DIR}/${KRAKEN_BASE}_${GROUP}_merged.sorted.bam" 2>>"$LOGFILE"
                
                # Calculer le taux de mapping pour les reads merged
                calculate_mapping_rate "${DAMAGE_DIR}/${KRAKEN_BASE}_${GROUP}_merged.sorted.bam" "$KRAKEN_BASE" "$GROUP" "merged"
                
                rm -f ${DAMAGE_DIR}/${KRAKEN_BASE}_${GROUP}_merged.sai ${DAMAGE_DIR}/${KRAKEN_BASE}_${GROUP}_merged.sam ${DAMAGE_DIR}/${KRAKEN_BASE}_${GROUP}_merged.bam 2>>"$LOGFILE"
                mapDamage -i "${DAMAGE_DIR}/${KRAKEN_BASE}_${GROUP}_merged.sorted.bam" -r "$REF_FASTA" --folder "${DAMAGE_DIR}/${KRAKEN_BASE}_${GROUP}_mapDamage_merged" 2>>"$LOGFILE"
            fi
        fi
    done
done

echo "Finished at $(date)" | tee -a "$LOGFILE"
echo "Mapping statistics saved in: $MAPPING_INFO" | tee -a "$LOGFILE"
