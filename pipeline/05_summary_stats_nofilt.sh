#!/usr/bin/env bash
#SBATCH -p short -N 1 -n 2 --mem 4gb --out logs/summary_stats_nofilt.log

module load prodigal
module unload miniconda2
module load miniconda3
# need python 3

GENOMES=genomes_nofilter
OUTPRED=results/prodigal_nofilter
OUTREPORT=results/genome_stats_nofilter
SEARCHDLV=search/prodigal_NCLDV_nofilter
SEARCHVOG=search/prodigal_NCVOG_nofilter
SEARCHPFAM=search/prodigal_pfamhits_nofilter
SEARCHRRNA=search/rRNA_nofilter
FUNGIRRNADB=db/rRNA/fungi.rRNA
SEARCHMT=search/MT_nofilter
FUNGIMTPEP=db/MT/MT_peps.fa
mkdir -p $OUTREPORT $OUTPRED $SEARCHRRNA $SEARCHMT
CPUS=$SLURM_CPUS_ON_NODE
if [ -z $CPUS ]; then
 CPUS=1
fi
CPU=$CPUS

INGENOME=$(ls $GENOMES/*.spades.fasta | sed -n ${N}p)
BASE=$(basename $INGENOME .spades.fasta)

for INGENOME in $(ls $GENOMES/*.spades.fasta)
do
  BASE=$(basename $INGENOME .spades.fasta)

  echo $BASE
  if [[ ! -f $OUTPRED/$BASE.prodigal.faa || $INGENOME -nt $OUTPRED/$BASE.prodigal.faa ]]; then
	   #prodigal -a $OUTPRED/$BASE.prodigal.faa -d $OUTPRED/$BASE.prodigal.cds -f gff -i $INGENOME -p meta -o $OUTPRED/$BASE.prodigal.gff
     echo "Expected to have already run prodigal and pre-screen Pfam? for $BASE ($OUTPRED/$BASE.prodigal.faa) or genome is newer than prodigal file"
     continue
   fi

   if [[ ! -s $SEARCHRRNA/$BASE.rRNA_search.tab ]]; then
     module load ncbi-blast/2.9.0+
     blastn -db $FUNGIRRNADB -num_threads $CPU -query $INGENOME -evalue 1e-8 -subject_besthit -outfmt 6 -max_target_seqs 5 -out $SEARCHRRNA/$BASE.rRNA_search.tab
     pigz -k $SEARCHRRNA/$BASE.rRNA_search.tab
     module unload ncbi-blast
   fi

   if [[ ! -f $SEARCHMT/$BASE.MT_search.tab || $INGENOME -nt $SEARCHMT/$BASE.MT_search.tab || $FUNGIMTPEP -nt $SEARCHMT/$BASE.MT_search.tab ]]; then
     module load fasta
     tfasty36 -E 1e-5 -m 8c -T $CPU $FUNGIMTPEP $INGENOME > $SEARCHMT/$BASE.MT_search.tab
     module unload fasta
   fi

   EXTRA=""
   if [ -s $SEARCHDLV/$BASE.hmmsearch.domtbl ]; then
     EXTRA="-sd $SEARCHDLV/$BASE.hmmsearch.domtbl"
   fi
   if [ -s $SEARCHVOG/$BASE.hmmsearch.domtbl ]; then
     EXTRA="$EXTRA -sv $SEARCHVOG/$BASE.hmmsearch.domtbl"
   fi
   if [ -s $SEARCHPFAM/$BASE.pfamscan.domtbl ]; then
     EXTRA="$EXTRA --pfam $SEARCHPFAM/$BASE.pfamscan.domtbl"
   fi
   python3 scripts/genome_stats_for_viral_ML.py --prodigal $OUTPRED/$BASE.prodigal.gff --outbase $OUTREPORT/$BASE \
     --ribosomal $SEARCHRRNA/$BASE.rRNA_search.tab --mito  $SEARCHMT/$BASE.MT_search.tab $EXTRA -i $INGENOME
done
