#!/usr/bin/env bash
#SBATCH -p short -N 1 -n 2 --mem 4gb --out logs/summary_stats_ref.%a.log

module load prodigal
module unload miniconda2
module load miniconda3
# need python 3

GENOMES=reference_genomes/DNA/
OUTPRED=results/prodigal_reference
OUTREPORT=results/genome_stats_ref
SEARCHDLV=search/genome_NCLDV
SEARCHVOG=search/genome_NCVOG
SEARCHRRNA=search/rRNA
FUNGIRRNADB=db/rRNA/fungi.rRNA
SEARCHMT=search/MT
FUNGIMTPEP=db/MT/MT_peps.fa
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

INGENOME=$(ls $GENOMES/*.nt.fasta | sed -n ${N}p)
BASE=$(basename $INGENOME .nt.fasta)

echo $BASE
if [ ! -f $INGENOME ]; then
  echo "cannot find a genome in $GENOMES for *.spades.fasta.gz"
  exit
fi
if [[ ! -f $OUTPRED/$BASE.prodigal.faa || $INGENOME -nt $OUTPRED/$BASE.prodigal.faa ]]; then
	prodigal -a $OUTPRED/$BASE.prodigal.faa -d $OUTPRED/$BASE.prodigal.cds -f gff -i $INGENOME -p meta -o $OUTPRED/$BASE.prodigal.gff
fi

if [[ ! -s $SEARCHRRNA/$BASE.rRNA_search.tab.gz ]]; then
  module load ncbi-blast/2.9.0+
  blastn -db $FUNGIRRNADB -num_threads $CPU -query $INGENOME -evalue 1e-3  -outfmt 6 -max_target_seqs 5 -out $SEARCHRRNA/$BASE.rRNA_search.tab
  pigz -k $SEARCHRRNA/$BASE.rRNA_search.tab
  module unload ncbi-blast
fi

if [[ ! -f $SEARCHMT/$BASE.MT_search.tab ]]; then
  module load fasta
  tfasty36 -E 1e-5 -m 8c -T $CPU $FUNGIMTPEP $INGENOME > $SEARCHMT/$BASE.MT_search.tab
  if [ ! -s $SEARCHMT/$BASE.MT_search.tab ]; then
    pigz -k $SEARCHMT/$BASE.MT_search.tab
  fi
  module unload fasta
fi

EXTRA=""
if [ -s $SEARCHDLV/$BASE.nt.TFASTX.tab ]; then
  EXTRA="$EXTRA -sd $SEARCHDLV/$BASE.nt.TFASTX.tab"
fi
if [ -s $SEARCHVOG/$BASE.nt.TFASTX.tab  ]; then
  EXTRA="$EXTRA -sv $SEARCHVOG/$BASE.nt.TFASTX.tab"
fi
python3 scripts/genome_stats_for_viral_ML.py --prodigal $OUTPRED/$BASE.prodigal.gff --outbase $OUTREPORT/$BASE --ribosomal $SEARCHRRNA/$BASE.rRNA_search.tab --mitochondria $SEARCHMT/$BASE.MT_search.tab $EXTRA -i $INGENOME
