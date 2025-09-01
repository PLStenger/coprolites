#!/bin/bash
#SBATCH --job-name=04_clumpify
#SBATCH --ntasks=1
#SBATCH -p smp
#SBATCH --mem=250G
#SBATCH --mail-user=pierrelouis.stenger@gmail.com
#SBATCH --mail-type=ALL 
#SBATCH --error="/home/plstenge/coprolites/00_scripts/04_clumpify.err"
#SBATCH --output="/home/plstenge/coprolites/00_scripts/04_clumpify.out"

# Dossiers
IN_DIR="/home/plstenge/coprolites/04_fastuniq"
OUT_DIR="/home/plstenge/coprolites/05_clumpify"

module load conda/4.12.0
source ~/.bashrc
conda activate bbduk

CLUMPIFY=/home/plstenge/bbmap/clumpify.sh

# Créer le dossier de sortie s'il n'existe pas
mkdir -p "$OUT_DIR"

# Boucle sur tous les fichiers R1
for R1 in "$IN_DIR"/*_R1.fastq; do
    # Déduire le nom du fichier R2 correspondant
    R2="${R1/_R1.fastq/_R2.fastq}"
    # Vérifier que le fichier R2 existe
    if [[ -f "$R2" ]]; then
        # Extraire le nom de base pour l'output
        base=$(basename "$R1" _R1.fastq)
        # Lancer clumpify pour la paire
        $CLUMPIFY \
            in="$R1" in2="$R2" \
            out="$OUT_DIR/${base}_clumpify_R1.fastq.gz" \
            out2="$OUT_DIR/${base}_clumpify_R2.fastq.gz" \
            dedupe=t
        echo "Traitement terminé pour $base"
    else
        echo "Fichier R2 manquant pour $R1, ignoré."
    fi
done
