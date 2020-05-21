#!/usr/bin/bash
#SBATCH -p short --ntasks 64 --mem 32gb --out logs/prodigal_pfam.%a.log

module load hmmer/3.3-mpi
module load db-pfam

CPUS=$SLURM_CPUS_ON_NODE
if [ -z $CPUS ]; then
 CPUS=1
fi
CPU=$CPUS

N=${SLURM_ARRAY_TASK_ID}

if [ -z $N ]; then
    N=$1
    if [ -z $N ]; then
        echo "Need an array id or cmdline val for the job"
        exit
    fi
fi
OUT=search/prodigal_pfamhits
GENOMES=genomes
OUTPRED=results/prodigal

FILE=$(ls $GENOMES/*.sorted.fasta | sed -n ${N}p)
BASE=$(basename $FILE .sorted.fasta)
TARGET=$OUTPRED/$BASE.prodigal.faa
PFAM=$PFAM_DB/Pfam-A.hmm
mkdir -p $OUT
OUTFILE=$OUT/${BASE}.pfamscan

if [[ ! -f $OUTFILE.done || $TARGET -nt $OUTFILE.done ]]; then
  srun hmmsearch --mpi --cut_ga --domtbl $OUTFILE.domtbl -o $OUTFILE.log $PFAM $TARGET
  touch $OUTFILE.done
fi
