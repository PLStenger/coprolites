#!/bin/bash
#SBATCH --job-name=09_bwa_n_test
#SBATCH --ntasks=1
#SBATCH -p smp
#SBATCH --mem=100G
#SBATCH --mail-user=pierrelouis.stenger@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --output=/home/plstenge/coprolites/00_scripts/09_bwa_n_test.out
#SBATCH --error=/home/plstenge/coprolites/00_scripts/09_bwa_n_test.err

module load conda/4.12.0
source ~/.bashrc
conda activate mapdamage_py39

# Paths
REF="/home/plstenge/genomes/Ovis_aries.ARS-UI_Ramb_v3.0.dna.toplevel.fa"
FASTQ="/home/plstenge/coprolites/06_fastp/clean_Aviti_cop408_dedup_clumpify_fastp_merged.fastq.gz"
OUTDIR="/home/plstenge/coprolites/08_damage/bwa_n_tests"
mkdir -p "$OUTDIR"

# bwa index once
bwa index "$REF"

# Different -n values to test
N_VALUES=("0" "0.02" "0.04" "0.06" "0.08" "1")

RESULTS="$OUTDIR/mapping_rates.txt"
echo -e "n\tTotalReads\tMappedReads\tMappingRate(%)" > "$RESULTS"

for n in "${N_VALUES[@]}"; do
    echo "Running bwa aln with -n $n ..."

    # Alignment
    bwa aln -n "$n" -l 24 -k 2 -q 20 -t 4 "$REF" "$FASTQ" > "$OUTDIR/aln_n${n}.sai"
    bwa samse "$REF" "$OUTDIR/aln_n${n}.sai" "$FASTQ" > "$OUTDIR/aln_n${n}.sam"

    # BAM conversion + sorting
    samtools view -bS "$OUTDIR/aln_n${n}.sam" | samtools sort -o "$OUTDIR/aln_n${n}.bam"
    samtools index "$OUTDIR/aln_n${n}.bam"
    rm "$OUTDIR/aln_n${n}.sam" "$OUTDIR/aln_n${n}.sai"

    # Get stats
    total=$(samtools view -c "$OUTDIR/aln_n${n}.bam")
    mapped=$(samtools view -c -F 4 "$OUTDIR/aln_n${n}.bam")
    rate=$(echo "scale=4; ($mapped/$total)*100" | bc)

    echo -e "${n}\t${total}\t${mapped}\t${rate}" >> "$RESULTS"
done

echo "All tests done. Results in $RESULTS"
