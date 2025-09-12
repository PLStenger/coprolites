#!/bin/bash
#SBATCH --job-name=99_adapterRemoval_whole_pipeline
#SBATCH --ntasks=1
#SBATCH -p smp
#SBATCH --mem=500G
#SBATCH --mail-user=pierrelouis.stenger@gmail.com
#SBATCH --mail-type=ALL 
#SBATCH --error="/home/plstenge/coprolites/00_scripts/99_adapterRemoval_whole_pipeline.err"
#SBATCH --output="/home/plstenge/coprolites/00_scripts/99_adapterRemoval_whole_pipeline.out"

# Fonction pour extraire le nom de base d'un échantillon
get_sample_basename() {
    local file="$1"
    # Extrait le nom sans les suffixes _R1/_R2 et extensions
    basename "$file" | sed -E 's/_R[12]\.(fastq|fq)(\.gz)?$//'
}

############################################################################################################
# 1) adapterremoval
############################################################################################################

WORKING_DIRECTORY="/home/plstenge/coprolites/01_raw_data"
OUTPUT="/home/plstenge/coprolites/99_AdapterRemoval_pipeline/02_cleaned_data_adapterremoval"
ADAPTER_FILE="/home/plstenge/seda_DNA_Corsican_wreck/99_softwares/adapters_adpateremoval.txt"

mkdir -p "$OUTPUT"

module load conda/4.12.0
source ~/.bashrc
conda activate adapterremoval

cd "$WORKING_DIRECTORY"

echo "=== Étape 1: AdapterRemoval - Démultiplexage et trimming ==="
for r1_file in *_R1.fastq.gz; do
    r2_file="${r1_file/_R1/_R2}"
    [[ ! -f "$r2_file" ]] && { echo "ERREUR: Fichier R2 manquant pour $r1_file" >&2; continue; }
    
    sample_name=$(get_sample_basename "$r1_file")
    echo "Traitement: $sample_name"
    
    AdapterRemoval \
        --adapter-list "$ADAPTER_FILE" \
        --file1 "$r1_file" \
        --file2 "$r2_file" \
        --basename "$OUTPUT/${sample_name}_step1" \
        --trimns \
        --trimqualities \
        --minlength 25 \
        --qualitymax 50 \
        --collapse \
        --gzip
done

conda deactivate

############################################################################################################
# 2) bbduk
############################################################################################################

INPUT_DIR="$OUTPUT"
OUTPUT="/home/plstenge/coprolites/99_AdapterRemoval_pipeline/03_bbduk"

mkdir -p "$OUTPUT"

module load conda/4.12.0
source ~/.bashrc
conda activate bbduk

PHIX=/home/plstenge/bbmap/resources/phix174_ill.ref.fa.gz
BBDUK=/home/plstenge/bbmap/bbduk.sh

cd "$INPUT_DIR"

echo "=== Étape 2: BBDuk - Suppression PhiX ==="
for r1_file in *_step1.pair1.fastq.gz; do
    r2_file="${r1_file/pair1/pair2}"
    [[ ! -f "$r2_file" ]] && { echo "ERREUR: Fichier pair2 manquant pour $r1_file" >&2; continue; }
    
    sample_name=$(echo "$r1_file" | sed 's/_step1\.pair1\.fastq\.gz$//')
    echo "Traitement: $sample_name"
    
    $BBDUK -Xmx4g \
        in1="$r1_file" \
        in2="$r2_file" \
        out1="$OUTPUT/${sample_name}_step2_R1.fastq.gz" \
        out2="$OUTPUT/${sample_name}_step2_R2.fastq.gz" \
        ref="$PHIX" \
        ktrim=rl \
        k=23 \
        mink=11 \
        hdist=1 \
        tpe \
        tbo \
        minlen=25 \
        qtrim=r \
        trimq=20 \
        stats="$OUTPUT/${sample_name}_step2_stats.txt"
done

conda deactivate

############################################################################################################
# 3) fastuniq
############################################################################################################

INPUT_DIR="$OUTPUT"
OUTPUT="/home/plstenge/coprolites/99_AdapterRemoval_pipeline/04_fastuniq"

mkdir -p "$OUTPUT"

module load conda/4.12.0
source ~/.bashrc
conda activate fastuniq

cd "$INPUT_DIR" || exit 1

TMP="/tmp/fastuniq_tmp"
mkdir -p "$TMP"

echo "=== Étape 3: FastUniq - Suppression des duplicatas ==="
for R1_gz in *_step2_R1.fastq.gz; do
    R2_gz="${R1_gz/_R1.fastq.gz/_R2.fastq.gz}"
    
    if [[ -f "$R2_gz" ]]; then
        sample_name=$(echo "$R1_gz" | sed 's/_step2_R1\.fastq\.gz$//')
        echo "Traitement: $sample_name"

        R1_tmp="${TMP}/${sample_name}_R1.fastq"
        R2_tmp="${TMP}/${sample_name}_R2.fastq"
        listfile="${TMP}/${sample_name}.list"

        zcat "$R1_gz" > "$R1_tmp"
        zcat "$R2_gz" > "$R2_tmp"

        echo -e "$R1_tmp\n$R2_tmp" > "$listfile"

        fastuniq -i "$listfile" -t q \
            -o "${OUTPUT}/${sample_name}_step3_R1.fastq" \
            -p "${OUTPUT}/${sample_name}_step3_R2.fastq"

        rm -f "$R1_tmp" "$R2_tmp" "$listfile"
    else
        echo "ATTENTION: fichier R2 manquant pour $R1_gz"
    fi
done

rm -rf "$TMP"
conda deactivate

############################################################################################################
# 4) clumpify
############################################################################################################

INPUT_DIR="$OUTPUT"
OUTPUT="/home/plstenge/coprolites/99_AdapterRemoval_pipeline/05_clumpify"

mkdir -p "$OUTPUT"

module load conda/4.12.0
source ~/.bashrc
conda activate bbduk

CLUMPIFY=/home/plstenge/bbmap/clumpify.sh

echo "=== Étape 4: Clumpify - Deduplication avancée ==="
for R1 in "$INPUT_DIR"/*_step3_R1.fastq; do
    R2="${R1/_R1.fastq/_R2.fastq}"
    
    if [[ -f "$R2" ]]; then
        sample_name=$(basename "$R1" _step3_R1.fastq)
        echo "Traitement: $sample_name"
        
        $CLUMPIFY \
            in="$R1" in2="$R2" \
            out="$OUTPUT/${sample_name}_step4_R1.fastq.gz" \
            out2="$OUTPUT/${sample_name}_step4_R2.fastq.gz" \
            dedupe=t
    else
        echo "Fichier R2 manquant pour $R1, ignoré."
    fi
done

conda deactivate

############################################################################################################
# 5) fastp
############################################################################################################

INPUT_DIR="$OUTPUT"
OUTPUT="/home/plstenge/coprolites/99_AdapterRemoval_pipeline/06_fastp"
LOG_DIR="/home/plstenge/coprolites/00_scripts/05_fastp_out/"

mkdir -p "$OUTPUT"
mkdir -p "$LOG_DIR"

module load conda/4.12.0
source ~/.bashrc
conda activate fastp

echo "=== Étape 5: Fastp - Nettoyage final et fusion ==="
for R1 in "$INPUT_DIR"/*_step4_R1.fastq.gz; do
    R2="${R1/_R1.fastq.gz/_R2.fastq.gz}"
    
    sample_name=$(basename "$R1" _step4_R1.fastq.gz)
    echo "Traitement: $sample_name"
    
    # Fichiers de sortie
    OUT_R1="$OUTPUT/${sample_name}_step5_R1.fastq.gz"
    OUT_R2="$OUTPUT/${sample_name}_step5_R2.fastq.gz"
    MERGED="$OUTPUT/${sample_name}_step5_merged.fastq.gz"
    HTML="$LOG_DIR/${sample_name}_fastp.html"
    JSON="$LOG_DIR/${sample_name}_fastp.json"
    
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
done

conda deactivate

############################################################################################################
# 6) kraken2
############################################################################################################

module load conda/4.12.0
source ~/.bashrc
conda activate kraken2

INPUT_DIR="$OUTPUT"
KRAKEN2_DB="/home/plstenge/k2_core_nt_20250609"
OUTPUT="/home/plstenge/coprolites/99_AdapterRemoval_pipeline/07_kraken2"
THREADS=36

mkdir -p "$OUTPUT"

echo "=== Étape 6: Kraken2 - Classification taxonomique ==="

# 1. Analyse des merged (single-end)
echo "Analyse Kraken2 des merged (single-end)"
for MERGED in "$INPUT_DIR"/*_step5_merged.fastq.gz; do
    sample_name=$(basename "$MERGED" _step5_merged.fastq.gz)
    echo "Traitement merged: $sample_name"
    
    OUT_KRAKEN="$OUTPUT/${sample_name}_merged.kraken"
    OUT_REPORT="$OUTPUT/${sample_name}_merged.report"

    kraken2 --conf 0.2 --db "$KRAKEN2_DB" --threads $THREADS \
        --output "$OUT_KRAKEN" --report "$OUT_REPORT" "$MERGED"
done

# 2. Analyse des unmerged (paired-end)
echo "Analyse Kraken2 des unmerged (paired-end)"
for R1 in "$INPUT_DIR"/*_step5_R1.fastq.gz; do
    sample_name=$(basename "$R1" _step5_R1.fastq.gz)
    R2="$INPUT_DIR/${sample_name}_step5_R2.fastq.gz"

    if [[ -f "$R2" ]]; then
        echo "Traitement unmerged: $sample_name"
        
        OUT_KRAKEN="$OUTPUT/${sample_name}_unmerged.kraken"
        OUT_REPORT="$OUTPUT/${sample_name}_unmerged.report"

        kraken2 --conf 0.2 --paired --db "$KRAKEN2_DB" --threads $THREADS \
            --output "$OUT_KRAKEN" --report "$OUT_REPORT" "$R1" "$R2"
    fi
done

conda deactivate

############################################################################################################
# 7) krona
############################################################################################################

module load conda/4.12.0
source ~/.bashrc
conda activate krona

INPUT_DIR="$OUTPUT"
OUTPUT="/home/plstenge/coprolites/99_AdapterRemoval_pipeline/08_krona"

mkdir -p "$OUTPUT"

echo "=== Étape 7: Krona - Visualisation taxonomique ==="

cd "$INPUT_DIR"
ktImportTaxonomy -t 5 -m 3 -o "$OUTPUT/multi-krona.html" "$INPUT_DIR"/*.report 

conda deactivate

############################################################################################################
# 8) MapDamage
############################################################################################################

module load conda/4.12.0
source ~/.bashrc
conda activate mapdamage_py39

KRAKEN_DIR="$INPUT_DIR"  # Répertoire des résultats Kraken2
FASTQ_DIR="/home/plstenge/coprolites/99_AdapterRemoval_pipeline/06_fastp"  # Répertoire des fichiers FASTQ traités
DAMAGE_BASE="/home/plstenge/coprolites/99_AdapterRemoval_pipeline/08_damage"
LOGFILE="/home/plstenge/coprolites/00_scripts/09_DAMAGE_$(date +%Y%m%d_%H%M%S).txt"
MAPPING_INFO="/home/plstenge/coprolites/99_AdapterRemoval_pipeline/08_damage/mapping_bwa_info.txt"

mkdir -p "$DAMAGE_BASE"

echo "=== Étape 8: MapDamage - Analyse des dommages ADN ancien ==="
echo "Script started at $(date)" | tee -a "$LOGFILE"

# Initialiser le fichier de mapping info avec les en-têtes
echo -e "Sample\tSpecies\tType\tTotal_Reads\tMapped_Reads\tMapping_Rate(%)" > "$MAPPING_INFO"

# Déclaration des génomes de référence
declare -A TAXONS=(
    ["Ovis_aries"]="9940:/home/plstenge/genomes/Ovis_aries.ARS-UI_Ramb_v3.0.dna.toplevel.fa"
    ["Capra_hircus"]="9925:/home/plstenge/genomes/Capra_hircus.ARS1.dna.toplevel.fa"
    ["Alnus_glutinosa"]="3517:/home/plstenge/genomes/Alnus_glutinosa_genome_assembly_dhAlnGlut1.fa"
    ["Corylus_avellana"]="13451:/home/plstenge/genomes/Corylus_avellana_CavTom2PMs_1_0.fasta"
)

# Fonction pour calculer le taux de mapping
calculate_mapping_rate() {
    local bam_file="$1"
    local sample_name="$2"
    local species="$3"
    local type="$4"
    
    if [[ -f "$bam_file" ]]; then
        local total_reads=$(samtools view -c "$bam_file")
        local mapped_reads=$(samtools view -c -F 4 "$bam_file")
        local mapping_rate=0
        if [[ $total_reads -gt 0 ]]; then
            mapping_rate=$(echo "scale=2; $mapped_reads * 100 / $total_reads" | bc)
        fi
        
        echo -e "${sample_name}\t${species}\t${type}\t${total_reads}\t${mapped_reads}\t${mapping_rate}" >> "$MAPPING_INFO"
        echo "Mapping stats for ${sample_name}_${species}_${type}: ${mapped_reads}/${total_reads} (${mapping_rate}%)" | tee -a "$LOGFILE"
    fi
}

shopt -s nullglob
for KRAKEN_FILE in "$KRAKEN_DIR"/*.kraken; do
    KRAKEN_BASE=$(basename "$KRAKEN_FILE" .kraken)
    echo -e "\n==== Processing file: $KRAKEN_FILE ====" | tee -a "$LOGFILE"
    
    # Extraire le nom de l'échantillon du fichier kraken
    if [[ "$KRAKEN_BASE" =~ (.+)_(merged|unmerged)$ ]]; then
        PREFIX="${BASH_REMATCH[1]}"
        TYPE="${BASH_REMATCH[2]}"
    else
        echo "Impossible d'extraire le préfixe de $KRAKEN_BASE" | tee -a "$LOGFILE"
        continue
    fi
    
    echo "Prefix: $PREFIX, Type: $TYPE" | tee -a "$LOGFILE"

    # Trouver les fichiers FASTQ correspondants
    if [[ "$TYPE" == "unmerged" ]]; then
        R1_FILES=("${FASTQ_DIR}/${PREFIX}_step5_R1.fastq"*)
        R2_FILES=("${FASTQ_DIR}/${PREFIX}_step5_R2.fastq"*)
        R1_FILE="${R1_FILES[0]:-}"
        R2_FILE="${R2_FILES[0]:-}"
        MERGED_FILE=""
    else
        MERGED_FILES=("${FASTQ_DIR}/${PREFIX}_step5_merged.fastq"*)
        MERGED_FILE="${MERGED_FILES[0]:-}"
        R1_FILE=""
        R2_FILE=""
    fi

    for GROUP in "${!TAXONS[@]}"; do
        TAX_ID="${TAXONS[$GROUP]%:*}"
        REF_FASTA="${TAXONS[$GROUP]#*:}"
        DAMAGE_DIR="${DAMAGE_BASE}/${GROUP}"
        mkdir -p "$DAMAGE_DIR"
        
        echo "Processing $GROUP (taxid: $TAX_ID) for $KRAKEN_BASE" | tee -a "$LOGFILE"
        
        OUT_R1="${DAMAGE_DIR}/${KRAKEN_BASE}_${GROUP}_R1.fastq"
        OUT_R2="${DAMAGE_DIR}/${KRAKEN_BASE}_${GROUP}_R2.fastq"
        OUT_MERGED="${DAMAGE_DIR}/${KRAKEN_BASE}_${GROUP}_merged.fastq"
        
        # Traitement des reads unmerged (paired-end)
        if [[ "$TYPE" == "unmerged" && -n "$R1_FILE" && -n "$R2_FILE" ]]; then
            echo "Extracting unmerged reads for $GROUP..." | tee -a "$LOGFILE"
            python3 /home/plstenge/KrakenTools/extract_kraken_reads.py \
                -k "$KRAKEN_FILE" -s "$R1_FILE" -s2 "$R2_FILE" -t "$TAX_ID" \
                -o "$OUT_R1" -o2 "$OUT_R2" --fastq-output 2>>"$LOGFILE"
                
            if [[ -f "$OUT_R1" && -f "$OUT_R2" && -s "$OUT_R1" && -s "$OUT_R2" ]]; then
                echo "Mapping unmerged reads for $GROUP..." | tee -a "$LOGFILE"
                bwa aln -n 0.08 -l 24 -k 2 -q 20 -t 4 "$REF_FASTA" "$OUT_R1" > "${DAMAGE_DIR}/${KRAKEN_BASE}_${GROUP}_R1.sai" 2>>"$LOGFILE"
                bwa aln -n 0.08 -l 24 -k 2 -q 20 -t 4 "$REF_FASTA" "$OUT_R2" > "${DAMAGE_DIR}/${KRAKEN_BASE}_${GROUP}_R2.sai" 2>>"$LOGFILE"
                bwa sampe "$REF_FASTA" "${DAMAGE_DIR}/${KRAKEN_BASE}_${GROUP}_R1.sai" "${DAMAGE_DIR}/${KRAKEN_BASE}_${GROUP}_R2.sai" "$OUT_R1" "$OUT_R2" > "${DAMAGE_DIR}/${KRAKEN_BASE}_${GROUP}.sam" 2>>"$LOGFILE"
                samtools view -bS "${DAMAGE_DIR}/${KRAKEN_BASE}_${GROUP}.sam" > "${DAMAGE_DIR}/${KRAKEN_BASE}_${GROUP}.bam" 2>>"$LOGFILE"
                samtools sort -o "${DAMAGE_DIR}/${KRAKEN_BASE}_${GROUP}.sorted.bam" "${DAMAGE_DIR}/${KRAKEN_BASE}_${GROUP}.bam" 2>>"$LOGFILE"
                samtools index "${DAMAGE_DIR}/${KRAKEN_BASE}_${GROUP}.sorted.bam" 2>>"$LOGFILE"
                
                calculate_mapping_rate "${DAMAGE_DIR}/${KRAKEN_BASE}_${GROUP}.sorted.bam" "$KRAKEN_BASE" "$GROUP" "unmerged"
                
                # Nettoyage des fichiers temporaires
                rm -f "${DAMAGE_DIR}/${KRAKEN_BASE}_${GROUP}_R1.sai" "${DAMAGE_DIR}/${KRAKEN_BASE}_${GROUP}_R2.sai" \
                      "${DAMAGE_DIR}/${KRAKEN_BASE}_${GROUP}.sam" "${DAMAGE_DIR}/${KRAKEN_BASE}_${GROUP}.bam" 2>>"$LOGFILE"
                
                # Analyse MapDamage
                mapDamage -i "${DAMAGE_DIR}/${KRAKEN_BASE}_${GROUP}.sorted.bam" \
                          -r "$REF_FASTA" \
                          --folder "${DAMAGE_DIR}/${KRAKEN_BASE}_${GROUP}_mapDamage_unmerged" 2>>"$LOGFILE"
            else
                echo "No reads extracted for $GROUP unmerged" | tee -a "$LOGFILE"
            fi
        fi
        
        # Traitement des reads merged (single-end)
        if [[ "$TYPE" == "merged" && -n "$MERGED_FILE" ]]; then
            echo "Extracting merged reads for $GROUP..." | tee -a "$LOGFILE"
            python3 /home/plstenge/KrakenTools/extract_kraken_reads.py \
                -k "$KRAKEN_FILE" -s "$MERGED_FILE" -t "$TAX_ID" -o "$OUT_MERGED" --fastq-output 2>>"$LOGFILE"
                
            if [[ -f "$OUT_MERGED" && -s "$OUT_MERGED" ]]; then
                echo "Mapping merged reads for $GROUP..." | tee -a "$LOGFILE"
                bwa aln -n 0.08 -l 24 -k 2 -q 20 -t 4 "$REF_FASTA" "$OUT_MERGED" > "${DAMAGE_DIR}/${KRAKEN_BASE}_${GROUP}_merged.sai" 2>>"$LOGFILE"
                bwa samse "$REF_FASTA" "${DAMAGE_DIR}/${KRAKEN_BASE}_${GROUP}_merged.sai" "$OUT_MERGED" > "${DAMAGE_DIR}/${KRAKEN_BASE}_${GROUP}_merged.sam" 2>>"$LOGFILE"
                samtools view -bS "${DAMAGE_DIR}/${KRAKEN_BASE}_${GROUP}_merged.sam" > "${DAMAGE_DIR}/${KRAKEN_BASE}_${GROUP}_merged.bam" 2>>"$LOGFILE"
                samtools sort -o "${DAMAGE_DIR}/${KRAKEN_BASE}_${GROUP}_merged.sorted.bam" "${DAMAGE_DIR}/${KRAKEN_BASE}_${GROUP}_merged.bam" 2>>"$LOGFILE"
                samtools index "${DAMAGE_DIR}/${KRAKEN_BASE}_${GROUP}_merged.sorted.bam" 2>>"$LOGFILE"
                
                calculate_mapping_rate "${DAMAGE_DIR}/${KRAKEN_BASE}_${GROUP}_merged.sorted.bam" "$KRAKEN_BASE" "$GROUP" "merged"
                
                # Nettoyage des fichiers temporaires
                rm -f "${DAMAGE_DIR}/${KRAKEN_BASE}_${GROUP}_merged.sai" "${DAMAGE_DIR}/${KRAKEN_BASE}_${GROUP}_merged.sam" \
                      "${DAMAGE_DIR}/${KRAKEN_BASE}_${GROUP}_merged.bam" 2>>"$LOGFILE"
                
                # Analyse MapDamage
                mapDamage -i "${DAMAGE_DIR}/${KRAKEN_BASE}_${GROUP}_merged.sorted.bam" \
                          -r "$REF_FASTA" \
                          --folder "${DAMAGE_DIR}/${KRAKEN_BASE}_${GROUP}_mapDamage_merged" 2>>"$LOGFILE"
            else
                echo "No reads extracted for $GROUP merged" | tee -a "$LOGFILE"
            fi
        fi
    done
done

echo "Finished at $(date)" | tee -a "$LOGFILE"
echo "Mapping statistics saved in: $MAPPING_INFO" | tee -a "$LOGFILE"

conda deactivate

echo "=== Pipeline terminé avec succès ==="



############################################################################################################
# 9) FastQC all steps
############################################################################################################

############################################################################################################
# 9) FastQC + MultiQC pour toutes les étapes
############################################################################################################

module load conda/4.12.0
source ~/.bashrc
conda activate fastqc  # Assurez-vous d'avoir un environnement avec fastqc et multiqc

BASE_DIR="/home/plstenge/coprolites"
QC_BASE_DIR="$BASE_DIR/99_AdapterRemoval_pipeline/09_quality_control"
THREADS=8

mkdir -p "$QC_BASE_DIR"

echo "=== Étape 9: FastQC + MultiQC - Contrôle qualité pour toutes les étapes ==="

# Définition des répertoires et étapes à analyser
declare -A STEPS_DIRS=(
    ["01_raw_data"]="$BASE_DIR/01_raw_data"
    ["02_adapterremoval"]="$BASE_DIR/99_AdapterRemoval_pipeline/02_cleaned_data_adapterremoval"
    ["03_bbduk"]="$BASE_DIR/99_AdapterRemoval_pipeline/03_bbduk"
    ["04_fastuniq"]="$BASE_DIR/99_AdapterRemoval_pipeline/04_fastuniq"
    ["05_clumpify"]="$BASE_DIR/99_AdapterRemoval_pipeline/05_clumpify"
    ["06_fastp"]="$BASE_DIR/99_AdapterRemoval_pipeline/06_fastp"
)

# Fonction pour exécuter FastQC sur un répertoire
run_fastqc_for_step() {
    local step_name="$1"
    local input_dir="$2"
    local output_dir="$QC_BASE_DIR/${step_name}_fastqc"
    
    echo "--- Analyse FastQC pour l'étape: $step_name ---"
    
    # Créer le répertoire de sortie
    mkdir -p "$output_dir"
    
    # Vérifier si le répertoire d'entrée existe
    if [[ ! -d "$input_dir" ]]; then
        echo "ATTENTION: Répertoire $input_dir n'existe pas, étape $step_name ignorée"
        return
    fi
    
    # Compter les fichiers FASTQ
    local fastq_count=$(find "$input_dir" -name "*.fastq*" -type f | wc -l)
    
    if [[ $fastq_count -eq 0 ]]; then
        echo "Aucun fichier FASTQ trouvé dans $input_dir"
        return
    fi
    
    echo "Traitement de $fastq_count fichiers FASTQ dans $input_dir"
    
    # Exécuter FastQC sur tous les fichiers FASTQ
    fastqc \
        --outdir "$output_dir" \
        --threads $THREADS \
        --quiet \
        "$input_dir"/*.fastq* 2>/dev/null
    
    echo "FastQC terminé pour $step_name ($fastq_count fichiers traités)"
}

# Exécuter FastQC pour chaque étape
for step_name in "${!STEPS_DIRS[@]}"; do
    run_fastqc_for_step "$step_name" "${STEPS_DIRS[$step_name]}"
done

# Créer un rapport MultiQC global
echo "--- Génération du rapport MultiQC global ---"
MULTIQC_OUTPUT="$QC_BASE_DIR/multiqc_report"
mkdir -p "$MULTIQC_OUTPUT"

# Exécuter MultiQC sur tous les résultats FastQC
multiqc \
    --outdir "$MULTIQC_OUTPUT" \
    --title "Pipeline AdapterRemoval - Rapport de qualité complet" \
    --comment "Analyse de qualité pour toutes les étapes du pipeline de traitement des données de coprolithes" \
    --dirs-depth 2 \
    --interactive \
    --export \
    "$QC_BASE_DIR"/*_fastqc/

echo "MultiQC terminé. Rapport disponible dans: $MULTIQC_OUTPUT/multiqc_report.html"

# Créer des rapports MultiQC séparés pour chaque étape (optionnel)
echo "--- Génération des rapports MultiQC par étape ---"
for step_name in "${!STEPS_DIRS[@]}"; do
    fastqc_dir="$QC_BASE_DIR/${step_name}_fastqc"
    
    if [[ -d "$fastqc_dir" ]] && [[ $(find "$fastqc_dir" -name "*.html" | wc -l) -gt 0 ]]; then
        multiqc_step_output="$QC_BASE_DIR/${step_name}_multiqc"
        mkdir -p "$multiqc_step_output"
        
        multiqc \
            --outdir "$multiqc_step_output" \
            --title "Étape $step_name - Rapport de qualité" \
            --filename "${step_name}_multiqc_report" \
            --interactive \
            "$fastqc_dir"
        
        echo "Rapport MultiQC créé pour l'étape $step_name"
    fi
done

# Créer un résumé des statistiques
echo "--- Génération du résumé des statistiques ---"
SUMMARY_FILE="$QC_BASE_DIR/pipeline_summary.txt"

cat > "$SUMMARY_FILE" << EOF
===========================================
RÉSUMÉ DU PIPELINE ADAPTERREMOVAL
===========================================
Date d'analyse: $(date)
Répertoire de base: $BASE_DIR

ÉTAPES ANALYSÉES:
EOF

for step_name in "${!STEPS_DIRS[@]}"; do
    input_dir="${STEPS_DIRS[$step_name]}"
    fastq_count=$(find "$input_dir" -name "*.fastq*" -type f 2>/dev/null | wc -l)
    echo "- $step_name: $fastq_count fichiers FASTQ" >> "$SUMMARY_FILE"
done

cat >> "$SUMMARY_FILE" << EOF

RAPPORTS GÉNÉRÉS:
- Rapport MultiQC global: $MULTIQC_OUTPUT/multiqc_report.html
- Rapports FastQC individuels: $QC_BASE_DIR/*_fastqc/
- Rapports MultiQC par étape: $QC_BASE_DIR/*_multiqc/

Pour visualiser les rapports:
1. Ouvrir le rapport MultiQC global dans un navigateur
2. Comparer les statistiques entre les étapes
3. Identifier les échantillons problématiques

===========================================
EOF

echo "Résumé sauvegardé dans: $SUMMARY_FILE"

# Afficher quelques statistiques rapides
echo ""
echo "=== RÉSUMÉ RAPIDE ==="
for step_name in "${!STEPS_DIRS[@]}"; do
    input_dir="${STEPS_DIRS[$step_name]}"
    if [[ -d "$input_dir" ]]; then
        fastq_count=$(find "$input_dir" -name "*.fastq*" -type f 2>/dev/null | wc -l)
        echo "$step_name: $fastq_count fichiers"
    fi
done

echo ""
echo "Contrôle qualité terminé avec succès !"
echo "Rapport principal: $MULTIQC_OUTPUT/multiqc_report.html"

conda deactivate
