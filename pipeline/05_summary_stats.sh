#!/usr/bin/env bash
#SBATCH -p short -N 1 -n 24 --mem 4gb --out logs/summary_stats.log

module load prodigal
module unload miniconda2
module load miniconda3
module unload perl
module load parallel
# need python 3

CPUS=$SLURM_CPUS_ON_NODE
if [ -z $CPUS ]; then
 CPUS=1
fi
CPU=$CPUS
export JOBCPU=2

run_sumstats() {
  FUNGIMTPEP=db/MT/MT_peps.fa
  FUNGIRRNADB=db/rRNA/fungi.rRNA
  echo "PRED=$OUTPRED"
  echo "REPORT=$OUTREPORT"
  echo "SEARCHDLV=$SEARCHDLV"
  INGENOME=$1
  EXT=$2
  BASE=$(basename $INGENOME $EXT)
  echo $BASE
  if [[ ! -f $OUTPRED/$BASE.prodigal.faa || $INGENOME -nt $OUTPRED/$BASE.prodigal.faa ]]; then
    #prodigal -a $OUTPRED/$BASE.prodigal.faa -d $OUTPRED/$BASE.prodigal.cds -f gff -i $INGENOME -p meta -o $OUTPRED/$BASE.prodigal.gff
    echo "Expected to have already run prodigal and pre-screen Pfam? for $BASE ($OUTPRED/$BASE.prodigal.faa) or genome is newer than prodigal file"
    continue
  fi

  if [[ ! -s $SEARCHRRNA/$BASE.rRNA_search.tab ]]; then
    module load ncbi-blast/2.9.0+
    blastn -db $FUNGIRRNADB -num_threads $JOBCPU -query $INGENOME -evalue 1e-8 -subject_besthit -outfmt 6 -max_target_seqs 5  -out $SEARCHRRNA/$BASE.rRNA_search.tab
    pigz -k $SEARCHRRNA/$BASE.rRNA_search.tab
    module unload ncbi-blast
  fi

  if [[ ! -f $SEARCHMT/$BASE.MT_search.tab || $INGENOME -nt $SEARCHMT/$BASE.MT_search.tab || $FUNGIMTPEP -nt $SEARCHMT/$BASE.MT_search.tab ]]; then
    module load fasta
    tfasty36 -E 1e-5 -m 8c -T $JOBCPU $FUNGIMTPEP $INGENOME > $SEARCHMT/$BASE.MT_search.tab
    module unload fasta
  fi
  # use the protein hmmsearch and GFF3

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

  #echo "searching for $SEARCHPFAM/$BASE.pfamscan.domtbl"
  #echo "EXTRA $EXTRA"
  python3 scripts/genome_stats_for_viral_ML.py --prodigal $OUTPRED/$BASE.prodigal.gff --outbase $OUTREPORT/$BASE \
    --ribosomal $SEARCHRRNA/$BASE.rRNA_search.tab -mt $SEARCHMT/$BASE.MT_search.tab $EXTRA -i $INGENOME
   echo in my_func $1
 }
export -f run_sumstats
# clean/sorted genomes
GENOMES=genomes
export OUTPRED=results/prodigal
export OUTREPORT=results/genome_stats
export SEARCHDLV=search/prodigal_NCLDV
export SEARCHVOG=search/prodigal_NCVOG
export SEARCHPFAM=search/prodigal_pfamhits
export SEARCHRRNA=search/rRNA
export SEARCHMT=search/MT

mkdir -p $OUTREPORT $OUTPRED $SEARCHRRNA $SEARCHMT

parallel --env run_sumstats --env JOBCPU --env SEARCHRRNA --env OUTPRED --env OUTREPORT --env SEARCHMT --env SEARCHDLV --env SEARCHVOG --env SEARCHPFAM -j $CPU run_sumstats {} .sorted.fasta ::: $(ls $GENOMES/*.sorted.fasta)


# reference genomes
GENOMES=reference_genomes/DNA/
export OUTPRED=results/prodigal_reference
export OUTREPORT=results/genome_stats_ref
export SEARCHDLV=search/prodigal_NCLDV
export SEARCHVOG=search/prodigal_NCVOG
export SEARCHRRNA=search/rRNA
export SEARCHMT=search/MT

mkdir -p $OUTREPORT $OUTPRED $SEARCHRRNA $SEARCHMT
parallel --env run_sumstats --env JOBCPU --env SEARCHRRNA --env OUTPRED --env OUTREPORT --env SEARCHMT --env SEARCHDLV --env SEARCHVOG --env SEARCHPFAM -j $CPU run_sumstats {} .nt.fasta ::: $(ls $GENOMES/*.nt.fasta)

# nofilter SPADes run
GENOMES=genomes_nofilter
export OUTPRED=results/prodigal_nofilter
export OUTREPORT=results/genome_stats_nofilter
export SEARCHDLV=search/prodigal_NCLDV_nofilter
export SEARCHVOG=search/prodigal_NCVOG_nofilter
export SEARCHPFAM=search/prodigal_pfamhits_nofilter
export SEARCHRRNA=search/rRNA_nofilter
export SEARCHMT=search/MT_nofilter

mkdir -p $OUTREPORT $OUTPRED $SEARCHRRNA $SEARCHMT
parallel --env run_sumstats --env JOBCPU --env SEARCHRRNA --env OUTPRED --env OUTREPORT --env SEARCHMT --env SEARCHDLV --env SEARCHVOG --env SEARCHPFAM -j $CPU run_sumstats {} .spades.fasta ::: $(ls $GENOMES/*.spades.fasta)
