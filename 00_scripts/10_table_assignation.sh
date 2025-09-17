#!/bin/bash
#SBATCH --job-name=10_table_assignation
#SBATCH --ntasks=1
#SBATCH -p smp
#SBATCH --mem=1000G
#SBATCH --mail-user=pierrelouis.stenger@gmail.com
#SBATCH --mail-type=ALL 
#SBATCH --error="/home/plstenge/coprolites/00_scripts/10_table_assignation.err"
#SBATCH --output="/home/plstenge/coprolites/00_scripts/10_table_assignation.out"

### cd /home/plstenge/coprolites/07_kraken2/
### 
### #git clone https://github.com/jenniferlu717/KrakenTools.git
### 
### python3 /home/plstenge/coprolites/07_kraken2/KrakenTools/kreport2mpa.py -r clean_Aviti_cop408_dedup_clumpify_merged.report  -o clean_Aviti_cop408_dedup_clumpify_merged.mpa
### python3 /home/plstenge/coprolites/07_kraken2/KrakenTools/kreport2mpa.py -r clean_Illumina_cop408_dedup_clumpify_merged.report  -o clean_Illumina_cop408_dedup_clumpify_merged.mpa
### python3 /home/plstenge/coprolites/07_kraken2/KrakenTools/kreport2mpa.py -r clean_Illumina_cop414_dedup_clumpify_merged.report  -o clean_Illumina_cop414_dedup_clumpify_merged.mpa
### python3 /home/plstenge/coprolites/07_kraken2/KrakenTools/kreport2mpa.py -r clean_Aviti_cop408_dedup_clumpify_unmerged.report  -o clean_Aviti_cop408_dedup_clumpify_unmerged.mpa
### python3 /home/plstenge/coprolites/07_kraken2/KrakenTools/kreport2mpa.py -r clean_Illumina_cop408_dedup_clumpify_unmerged.report  -o clean_Illumina_cop408_dedup_clumpify_unmerged.mpa
### python3 /home/plstenge/coprolites/07_kraken2/KrakenTools/kreport2mpa.py -r clean_Illumina_cop414_dedup_clumpify_unmerged.report  -o clean_Illumina_cop414_dedup_clumpify_unmerged.mpa
### python3 /home/plstenge/coprolites/07_kraken2/KrakenTools/kreport2mpa.py -r clean_Aviti_cop412_dedup_clumpify_merged.report  -o clean_Aviti_cop412_dedup_clumpify_merged.mpa
### python3 /home/plstenge/coprolites/07_kraken2/KrakenTools/kreport2mpa.py -r clean_Illumina_cop410_dedup_clumpify_merged.report  -o clean_Illumina_cop410_dedup_clumpify_merged.mpa
### python3 /home/plstenge/coprolites/07_kraken2/KrakenTools/kreport2mpa.py -r clean_Illumina_cop417_dedup_clumpify_merged.report  -o clean_Illumina_cop417_dedup_clumpify_merged.mpa
### python3 /home/plstenge/coprolites/07_kraken2/KrakenTools/kreport2mpa.py -r clean_Aviti_cop412_dedup_clumpify_unmerged.report  -o clean_Aviti_cop412_dedup_clumpify_unmerged.mpa
### python3 /home/plstenge/coprolites/07_kraken2/KrakenTools/kreport2mpa.py -r clean_Illumina_cop410_dedup_clumpify_unmerged.report  -o clean_Illumina_cop410_dedup_clumpify_unmerged.mpa
### python3 /home/plstenge/coprolites/07_kraken2/KrakenTools/kreport2mpa.py -r clean_Illumina_cop417_dedup_clumpify_unmerged.report  -o clean_Illumina_cop417_dedup_clumpify_unmerged.mpa
### python3 /home/plstenge/coprolites/07_kraken2/KrakenTools/kreport2mpa.py -r clean_Aviti_cop414_dedup_clumpify_merged.report  -o clean_Aviti_cop414_dedup_clumpify_merged.mpa
### python3 /home/plstenge/coprolites/07_kraken2/KrakenTools/kreport2mpa.py -r clean_Illumina_cop412_dedup_clumpify_merged.report  -o clean_Illumina_cop412_dedup_clumpify_merged.mpa
### python3 /home/plstenge/coprolites/07_kraken2/KrakenTools/kreport2mpa.py -r clean_Illumina_NTC_cop_dedup_clumpify_merged.report  -o clean_Illumina_NTC_cop_dedup_clumpify_merged.mpa
### python3 /home/plstenge/coprolites/07_kraken2/KrakenTools/kreport2mpa.py -r clean_Aviti_cop414_dedup_clumpify_unmerged.report  -o clean_Aviti_cop414_dedup_clumpify_unmerged.mpa
### python3 /home/plstenge/coprolites/07_kraken2/KrakenTools/kreport2mpa.py -r clean_Illumina_cop412_dedup_clumpify_unmerged.report  -o clean_Illumina_cop412_dedup_clumpify_unmerged.mpa
### python3 /home/plstenge/coprolites/07_kraken2/KrakenTools/kreport2mpa.py -r clean_Illumina_NTC_cop_dedup_clumpify_unmerged.report  -o clean_Illumina_NTC_cop_dedup_clumpify_unmerged.mpa
### 
### python3 /home/plstenge/coprolites/07_kraken2/KrakenTools/combine_mpa.py -i clean_Aviti_cop408_dedup_clumpify_merged.mpa clean_Aviti_cop408_dedup_clumpify_unmerged.mpa clean_Aviti_cop412_dedup_clumpify_merged.mpa clean_Aviti_cop412_dedup_clumpify_unmerged.mpa clean_Aviti_cop414_dedup_clumpify_merged.mpa clean_Aviti_cop414_dedup_clumpify_unmerged.mpa clean_Illumina_cop408_dedup_clumpify_merged.mpa clean_Illumina_cop408_dedup_clumpify_unmerged.mpa clean_Illumina_cop410_dedup_clumpify_merged.mpa clean_Illumina_cop410_dedup_clumpify_unmerged.mpa clean_Illumina_cop412_dedup_clumpify_merged.mpa clean_Illumina_cop412_dedup_clumpify_unmerged.mpa clean_Illumina_cop414_dedup_clumpify_merged.mpa clean_Illumina_cop414_dedup_clumpify_unmerged.mpa clean_Illumina_cop417_dedup_clumpify_merged.mpa clean_Illumina_cop417_dedup_clumpify_unmerged.mpa clean_Illumina_NTC_cop_dedup_clumpify_merged.mpa clean_Illumina_NTC_cop_dedup_clumpify_unmerged.mpa -o combined_mpa.tsv

#!/bin/bash
cd /home/plstenge/coprolites/99_AdapterRemoval_pipeline/07_kraken2

for file in *.report; do
  base=$(basename "$file" .report)
  python3 /home/plstenge/coprolites/07_kraken2/KrakenTools/kreport2mpa.py -r "$file" -o "${base}.mpa"
done

for file in *.mpa; do
  name=$(basename "$file" .mpa)
  sed -i "1i #SampleName\t${name}" "$file"
done

python3 /home/plstenge/coprolites/07_kraken2/KrakenTools/combine_mpa.py -i *.mpa -o combined_mpa.tsv

