#!/bin/bash
#SBATCH --job-name=99_goat_genome_kraken
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=36
#SBATCH -p smp
#SBATCH --mem=500G
#SBATCH --mail-user=pierrelouis.stenger@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --error="/home/plstenge/coprolites/00_scripts/goat_kraken.err"
#SBATCH --output="/home/plstenge/coprolites/00_scripts/goat_kraken.out"

# Activation environnement
module load conda/4.12.0
source ~/.bashrc
conda activate kraken2

# Répertoires
GENOME="/home/plstenge/genomes/Capra_hircus.ARS1.dna.toplevel.fa"
WORKDIR="/home/plstenge/coprolites/07_goat_kraken"
KRAKEN2_DB="/home/plstenge/k2_core_nt_20250609"
THREADS=36

mkdir -p "$WORKDIR"
cd "$WORKDIR"

########################################
# 1. Fractionner le génome en fragments
########################################
# On coupe en fenêtres de 150 bases, chevauchement 0
echo "Découpage du génome en faux reads..."
module load seqkit
seqkit sliding --window 150 --step 150 "$GENOME" > goat_genome_reads.fa

########################################
# 2. Conversion éventuelle en fastq
########################################
# Kraken2 accepte fasta ou fastq. Donc fastq n’est pas obligatoire.
# Mais si besoin d'un fastq fictif avec une qualité par défaut :
# awk 'NR%2==1 {print "@"substr($0,2)} NR%2==0 {print $0"\n+\n"gensub(/./,"I","g",$0)}' goat_genome_reads.fa > goat_genome_reads.fq

########################################
# 3. Lancer Kraken2 sur ces fragments
########################################
echo "Lancement Kraken2..."
kraken2 --db "$KRAKEN2_DB" \
        --threads $THREADS \
        --conf 0.2 \
        --output goat_genome.kraken \
        --report goat_genome.report \
        goat_genome_reads.fa

echo "Analyse terminée. Résultats : goat_genome.report"
