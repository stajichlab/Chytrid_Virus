#!/usr/bin/env bash
#SBATCH -p short -N 1 -n 2 --mem 4gb --out logs/summary_stats.%a.log

module load prodigal
module unload miniconda2
module load miniconda3
# need python 3

GENOMES=genomes_nofilter
OUTPRED=results/prodigal_nofilter
OUTREPORT=results/genome_stats_nofilter

mkdir -p $OUTREPORT $OUTPRED
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

INGENOME=$(ls $GENOMES/*.spades.fasta | sed -n ${N}p)
BASE=$(basename $INGENOME .spades.fasta)

echo $BASE
if [ ! -f $INGENOME ]; then
  echo "cannot find a genome in $GENOMES for *.spades.fasta.gz"
  exit
fi
if [[ ! -f $OUTPRED/$BASE.prodigal.faa || $INGENOME -nt $OUTPRED/$BASE.prodigal.faa ]]; then
	prodigal -a $OUTPRED/$BASE.prodigal.faa -d $OUTPRED/$BASE.prodigal.cds -f gff -i $INGENOME -p meta -o $OUTPRED/$BASE.prodigal.gff
fi

python3 scripts/genome_stats_for_viral_ML.py --prodigal $OUTPRED/$BASE.prodigal.gff --outbase $OUTREPORT/$BASE -i $INGENOME
