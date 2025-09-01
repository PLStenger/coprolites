#!/bin/bash
#SBATCH --job-name=06_kraken_db
##SBATCH --time=96:00:00       
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8     
#SBATCH -p smp
#SBATCH --mem=1000G            
#SBATCH --mail-user=pierrelouis.stenger@gmail.com
#SBATCH --mail-type=ALL 
#SBATCH --error="/home/plstenge/coprolites/00_scripts/06_kraken_db.err"
#SBATCH --output="/home/plstenge/coprolites/00_scripts/06_kraken_db.out"

module load conda/4.12.0
source ~/.bashrc
conda activate kraken2


# Aller sur https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/
# et télécharger:
# nucl_gb.accession2taxid.gz
# nucl_wgs.accession2taxid.gz
# pdb.accession2taxid.gz
# dead_nucl.accession2taxid.gz
# ex: wget https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/dead_nucl.accession2taxid.gz

# gunzip nucl_gb.accession2taxid.gz
# gunzip nucl_wgs.accession2taxid.gz
# gunzip pdb.accession2taxid.gz
# gunzip dead_nucl.accession2taxid.gz

# mv nucl_gb.accession2taxid /home/plstenge/nt_kraken2_db/taxonomy/
# mv nucl_wgs.accession2taxid /home/plstenge/nt_kraken2_db/taxonomy/
# mv pdb.accession2taxid /home/plstenge/nt_kraken2_db/taxonomy/
# mv dead_nucl.accession2taxid /home/plstenge/nt_kraken2_db/taxonomy/

# Créer dossier local
mkdir -p /home/plstenge/nt_kraken2_db_big
mkdir -p /home/plstenge/nt_kraken2_db_parts

# Diviser le FASTA en 10 fichiers
seqkit split -p 10 /storage/biodatabanks/ncbi/NT/ncbi_blast_nt_2024-8-24/fasta/All/all.fasta -O /home/plstenge/nt_kraken2_db_parts/

# Ajouter chaque morceau à la base
for f in /home/plstenge/nt_kraken2_db_parts/*.fasta; do
  kraken2-build --add-to-library "$f" --db /home/plstenge/nt_kraken2_db_big
done

# Télécharger taxonomie
kraken2-build --download-taxonomy --db /home/plstenge/nt_kraken2_db_big

# Construire la base
kraken2-build --build --db /home/plstenge/nt_kraken2_db_big --threads 36


#kraken2-build --add-to-library /storage/biodatabanks/ncbi/NT/ncbi_blast_nt_2024-8-24/fasta/All/all.fasta --db /home/plstenge/nt_kraken2_db

#echo "add-to-library ok"

#kraken2-build --download-taxonomy --db /home/plstenge/nt_kraken2_db

#echo "download-taxonomy ok"

#kraken2-build --build --db /home/plstenge/nt_kraken2_db --threads 36
#kraken2-build --build --db /home/plstenge/nt_kraken2_db --threads 36 > build.log 2> build.err

#echo "base de donnée construite"
