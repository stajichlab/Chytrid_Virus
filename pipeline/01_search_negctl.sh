#!/usr/bin/bash
#SBATCH -p intel -N 1 -n 8 --mem 32gb --out logs/search_netctl.%a.log

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
GENOMES=negativectl_genomes
QUERY=db/NCVOG/NCVOG.consensus.fasta
OUT=search/negctl_genome
mkdir -p $OUT
TARGET=$(ls $GENOMES/*.nt.fasta | sed -n ${N}p)
OUTFILE=$OUT/$(basename $TARGET .nt.fasta).TFASTX.tab

if [[ ! -f $OUTFILE.done || $TARGET -nt $OUTFILE.done ]]; then
  tfasty36 -m 8c -E 1e-8 -T $CPU $QUERY $TARGET > $OUTFILE
  touch $OUTFILE.done
fi
QUERY=db/Guglielmini/NCLDV.consensus.fasta
OUT=search/negctl_genome_NCLDV
mkdir -p $OUT
OUTFILE=$OUT/$(basename $TARGET .fasta).TFASTX.tab
if [[ ! -f $OUTFILE.done || $TARGET -nt $OUTFILE.done ]]; then
	tfasty36 -m 8c -E 1e-10 -T $CPU $QUERY $TARGET > $OUTFILE
	touch $OUTFILE.done
fi

