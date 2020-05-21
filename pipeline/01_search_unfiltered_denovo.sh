#!/usr/bin/bash
#SBATCH -p intel -N 1 -n 8 --mem 32gb --out logs/search_unfiltered.%a.log
module load fasta

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
QUERY=db/NCVOG/NCVOG.consensus.fasta
OUT=search/genome_NCVOG
OUTPRED=results/prodigal_nofilter
mkdir -p $OUT
TARGET=$(ls $GENOMES/*.spades.fasta | sed -n ${N}p)


if [[ ! -f $OUTPRED/$BASE.prodigal.faa || $INGENOME -nt $OUTPRED/$BASE.prodigal.faa ]]; then
	prodigal -a $OUTPRED/$BASE.prodigal.faa -d $OUTPRED/$BASE.prodigal.cds -f gff -i $INGENOME -p meta -o $OUTPRED/$BASE.prodigal.gff
fi

if [[ ! -f $OUTFILE.done || $TARGET -nt $OUTFILE.done ]]; then
  tfasty36 -m 8c -E 1e-10 -T $CPU $QUERY $TARGET > $OUTFILE
  touch $OUTFILE.done
fi

QUERY=db/Guglielmini/NCLDV.consensus.fasta
OUT=search/genome_NCLDV
mkdir -p $OUT
OUTFILE=$OUT/$(basename $TARGET .fasta).TFASTX.tab
if [[ ! -f $OUTFILE.done || $TARGET -nt $OUTFILE.done ]]; then
	tfasty36 -m 8c -E 1e-10 -T $CPU $QUERY $TARGET > $OUTFILE
	touch $OUTFILE.done
fi
