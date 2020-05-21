#!/usr/bin/bash
#SBATCH -p short --ntasks 32 --mem 32gb --out logs/prodigal_pfam.log

module load hmmer/3.3-mpi
module load db-pfam

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
  echo
  hmmsearch --mpi --cut_ga --domtbl $OUTFILE.domtbl -o $OUTFILE.log $PFAM $TARGET
  touch $OUTFILE.done
fi
