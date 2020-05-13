#!/usr/bin/bash
#SBATCH -p short

if [ ! -d spades_nofilter ]; then
	# replace with a osf.io download in future?
	gdrive download -r 1Bv1MkQzki9GJLpE8iziCSaUuhD5hk1_s
	mv spades_noreadfilter_for_Virus spades_nofilter
	pushd genomes_nofilter
	for file in *.gz; do m=$(echo $file | perl -p -e 's/^\S+_((JEL|ARG|KLL_TL\S+|CCIBt|Burma_1F|California_|PLAUS|WJD)\d*)\.spades/$1.spades/'); mv $file $m; done
	popd
	pigz -d genomes_nofilter/*.gz
fi
pushd db/NCVOG

# next make the VOG files
# or replace with
#bash make_NCVOG_alns.sh
sbatch make_NCVOG_alns.sh

popd
mkdir -p db/Guglielmini
pushd db/Guglielmini
curl -O "https://zenodo.org/record/3368642/files/Additional%20data.zip"
unzip Additional%20data.zip

mkdir -p db/rRNA
pushd db/rRNA
module load ncbi-blast/2.9.0+
if [ ! -s fungi.rRNA.fasta ]; then
	  for FILE in ITS 28SrRNA 18SrRNA; do
			curl -O -L -C - "ftp://ftp.ncbi.nlm.nih.gov/refseq/TargetedLoci/Fungi/fungi.$FILE.fna.gz"
		done
		pigz -dc fungi.*.fna.gz > fungi.rRNA.fasta
		makeblastdb -in fungi.rRNA.fasta -title fungi.rRNA -out fungi.rRNA -dbtype nucl -parse_seqids
		rm -f *.gz
fi
popd

mkdir -p db/MT
pushd db/MT
if [ ! -s MT_peps.fa ]; then
	for N in P80440 Q0H8Y4 P15563 P03880 P50363 P50365 P50368 P80439 P80441 Q37395 ; do
		curl "https://www.uniprot.org/uniprot/P80440.fasta" >> MT_peps.fa
	done
fi
popd
mkdir -p gff3
cat lib/ncbi_download.txt  | while read GFF URL
do
	curl $URL | zgrep -P "\tgene\t" | perl -p -e 's/ID=gene-([^;]+);.+/ID=$1/' | pigz -c > gff3/$GFF.gz
done
