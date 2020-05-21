#!/usr/bin/bash
#SBATCH -p short -N 1 -n 4 --mem 2gb --out logs/prodigal_pep_search_unfil.%a.log

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

GENOMES=genomes_nofilter
OUTPRED=results/prodigal_nofilter
TARGET=$(ls $GENOMES/*.spades.fasta | sed -n ${N}p)
BASE=$(basename $TARGET .spades.fasta)
TARGET=$OUTPRED/$BASE.prodigal.faa
if [[ ! -f $TARGET || $INGENOME -nt $TARGET ]]; then
  if [ -f $TARGET.gz ]; then
    pigz -k $OUTPRED/$BASE.prodigal.*.gz
  else
	   prodigal -a $TARGET -d $OUTPRED/$BASE.prodigal.cds -f gff -i $INGENOME -p meta -o $OUTPRED/$BASE.prodigal.gff
   fi
fi

OUT=search/prodigal_NCLDV_nofilter
QUERY=db/Guglielmini/NCLDV.hmmb
mkdir -p $OUT
OUTFILE=$OUT/${BASE}.hmmsearch

if [[ ! -f $OUTFILE.done || $TARGET -nt $OUTFILE.done ]]; then
        hmmsearch -E 1e-8 --cpu $CPU --domtblout $OUTFILE.domtbl $QUERY $TARGET > $OUTFILE.log
        touch $OUTFILE.done
fi

QUERY=db/NCVOG/NCVOG.hmmb
OUT=search/prodigal_NCVOG_nofilter
mkdir -p $OUT

OUTFILE=$OUT/${BASE}.hmmsearch

if [[ ! -f $OUTFILE.done || $TARGET -nt $OUTFILE.done ]]; then
	hmmsearch -E 1e-8 --cpu $CPU --domtblout $OUTFILE.domtbl $QUERY $TARGET > $OUTFILE.log
  touch $OUTFILE.done
fi
