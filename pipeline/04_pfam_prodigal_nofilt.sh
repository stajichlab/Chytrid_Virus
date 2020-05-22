#!/usr/bin/bash
#SBATCH  -p intel --ntasks 64 --mem 24gb --out logs/prodigal_pfam_nofilt.%a.log --time 12:00:00

module load hmmer/3.3-mpi
module load db-pfam

N=${SLURM_ARRAY_TASK_ID}

if [ -z $N ]; then
    N=$1
    if [ -z $N ]; then
        echo "Need an array id or cmdline val for the job"
        exit
    fi
fi

OUT=search/prodigal_pfamhits_nofilter
GENOMES=genomes_nofilter
OUTPRED=results/prodigal_nofilter

FILE=$(ls $GENOMES/*.spades.fasta | sed -n ${N}p)
BASE=$(basename $FILE .spades.fasta)
TARGET=$OUTPRED/$BASE.prodigal.faa
PFAM=$PFAM_DB/Pfam-A.hmm
mkdir -p $OUT
OUTFILE=$OUT/${BASE}.pfamscan

if [[ ! -f $OUTFILE.done || $TARGET -nt $OUTFILE.done ]]; then
  time srun hmmsearch --mpi --cut_ga --domtbl $OUTFILE.domtbl -o $OUTFILE.log $PFAM $TARGET
fi
