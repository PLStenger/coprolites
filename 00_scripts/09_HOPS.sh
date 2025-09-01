#!/bin/bash
#SBATCH --job-name=09_HOPS
#SBATCH --ntasks=1
#SBATCH -p smp
#SBATCH --mem=1000G
#SBATCH --mail-user=pierrelouis.stenger@gmail.com
#SBATCH --mail-type=ALL 
#SBATCH --error="/home/plstenge/coprolites/00_scripts/09_HOPS.err"
#SBATCH --output="/home/plstenge/coprolites/00_scripts/09_HOPS.out"

module load conda/4.12.0
source ~/.bashrc
conda activate hops

# Modified script originaly from Buchwald et al 2024 - Plankton community changes during the last 124 000 years

# set variables (requires modification)
#================================================
WORKDIR=/home/plstenge/coprolites/06_fastp
INPUT=/home/plstenge/coprolites/06_fastp/clean_1131_sed8_rep3_dedup_clumpify_fastp_merged.fastq.gz
OUTPUT=/home/plstenge/coprolites/09_damages
CONFIG=/home/plstenge/coprolites/00_scripts/09_HOPS_configfile_nt.txt
INDEXDB=/path/to/database/nt-hops-20-11-03_step8/
TAXFILE=path/to/working/directory/taxalist_2022_09_20.txt
NCBIRESC=/path/to/database/ncbi
MEM=1400

mkdir -p $OUTPUT

# preparing the working environment
#================================================
cd $WORKDIR
module load bio/hops/0.34

# create Hops config file
#================================================
echo "preProcess=0 
alignment=1 
dupRemOff=0 
filter=def_anc 
index=${INDEXDB}
pathToList=${TAXFILE}
resources=${NCBIRESC}

useSlurm = 0
threadsMalt=${SLURM_CPUS_PER_TASK}
maxMemoryMalt=${MEM}

threadsMaltEx=${SLURM_CPUS_PER_TASK}
maxMemoryMaltEx=${MEM}

threadsPost=${SLURM_CPUS_PER_TASK}
maxMemoryPost=256" > ${CONFIG}

# tasks to be performed
#================================================

srun hops -Xmx${MEM}G -input ${INPUT} -output ${OUTPUT} -m full -c ${CONFIG} 


