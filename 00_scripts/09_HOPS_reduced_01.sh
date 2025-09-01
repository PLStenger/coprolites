#!/bin/bash
#SBATCH --job-name=09_HOPS_reduced_01
#SBATCH --ntasks=1
#SBATCH -p smp
#SBATCH --mem=1000G
#SBATCH --mail-user=pierrelouis.stenger@gmail.com
#SBATCH --mail-type=ALL 
#SBATCH --error="/home/plstenge/coprolites/00_scripts/09_HOPS_reduced_01.err"
#SBATCH --output="/home/plstenge/coprolites/00_scripts/09_HOPS_reduced_01.out"

#set -euo pipefail

KR_DIR="/home/plstenge/coprolites/07_kraken2"
OUT_DIR="/home/plstenge/custom_nt_hopped"
ACC2TAXID="/home/plstenge/custom_nt_hopped/nucl_gb.accession2taxid"

module load conda/4.12.0
source ~/.bashrc
conda activate hops

mkdir -p $OUT_DIR
cd $OUT_DIR

##########################################################################################
# À faire une seule fois
##########################################################################################
#mkdir /home/plstenge/custom_nt_hopped
#cd /home/plstenge/custom_nt_hopped
#
## Récupérer le mapping accession2taxid 
#wget https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz
#gunzip nucl_gb.accession2taxid.gz
#
## Extraire les accessions de all.fasta (NT - NCBI)
#grep '^>' /storage/biodatabanks/ncbi/NT/current/fasta/All/all.fasta | sed 's/^>\([^.]*\).*/\1/' | sort -u > all_accessions.txt
#
##  Filtrer le mapping accession2taxid
#awk 'NR==FNR{a[$1];next} ($1 in a)' all_accessions.txt ${ACC2TAXID} > filtered_accession2taxid.txt
##########################################################################################

### 2. Extraire les taxons détectés par Kraken2 (tous les rangs)
# Concatenation de tous les fichiers .report
# Extraction des taxIDs (script Python rapide)
#cat ${KR_DIR}/*.report > all_kraken_reports.txt

#awk -F"\t" '{if($5 ~ /^[0-9]+$/) print $5}' all_kraken_reports.txt | sort -u > detected_taxids.txt

# Filtrer le mapping accession2taxid pour ne garder que nos taxons
# Convertir le mapping en {taxid: [accession1, accession2, ...]}
#awk -F"\t" 'NR==FNR{a[$1]; next} ($3 in a)' detected_taxids.txt filtered_accession2taxid.txt > final_accession2taxid.txt

# Extraire les lignes de all.fasta qui correspondent à nos accessions filtrées
# Liste des accessions à garder
#cut -f1 final_accession2taxid.txt | sort -u > selected_accessions.txt

# Créer un fasta filtré (sélectionne les séquences correspondantes)
python3 -c "
import sys
accs = set()
with open('selected_accessions.txt') as f:
    for line in f:
        accs.add(line.strip())
selecting = False
print('Processing FASTA...', file=sys.stderr)
with open('/storage/biodatabanks/ncbi/NT/current/fasta/All/all.fasta') as f:
    for line in f:
        if line.startswith('>'):
            acc = line[1:].split('.')[0].split()[0]
            if acc in accs:
                selecting = True
                print(line, end='')
            else:
                selecting = False
        elif selecting:
            print(line, end='')
" > custom_nt.fasta
