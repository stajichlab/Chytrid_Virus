#!/usr/bin/bash
#SBATCH -p short -N 1 -n 2 --mem 4gb --out logs/pepsearch.%a.log

module load hmmer

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
PEP=proteins
QUERY=db/NCVOG/NCVOG.hmmb
OUT=search/protein_NCVOG
mkdir -p $OUT
TARGET=$(ls $PEP/*.aa.fasta | sed -n ${N}p)
OUTFILE=$OUT/$(basename $TARGET .aa.fasta).hmmscan

if [[ ! -f $OUTFILE.done || $TARGET -nt $OUTFILE.done ]]; then
	hmmscan -E 1e-8 --cpu $CPU --domtblout $OUTFILE.domtbl $QUERY $TARGET > $OUTFILE.log
  touch $OUTFILE.done
fi

OUT=search/protein_NCLDV
QUERY=db/Guglielmini/NCLDV.hmmb
mkdir -p $OUT
OUTFILE=$OUT/$(basename $TARGET .aa.fasta).hmmscan
if [[ ! -f $OUTFILE.done || $TARGET -nt $OUTFILE.done ]]; then
        hmmscan -E 1e-8 --cpu $CPU --domtblout $OUTFILE.domtbl $QUERY $TARGET > $OUTFILE.log
  touch $OUTFILE.done
fi
