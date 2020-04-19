#!/usr/bin/bash
#SBATCH -p short -N 1 -n 24 --mem 32gb 

module unload perl
module load parallel
module load ncbi-blast/2.8.0+
module load muscle
module load hmmer/3
FORCE=0
CPUS=$SLURM_CPUS_ON_NODE
if [ -z $CPUS ]; then
 CPUS=1
fi
CPU=$CPUS

if [ ! -f NCVOG.csv ]; then
	curl -O ftp://ftp.ncbi.nih.gov/pub/wolf/COGs/NCVOG/NCVOG.csv
fi
if [ ! -f NCVOG.fa ]; then
	curl -O ftp://ftp.ncbi.nih.gov/pub/wolf/COGs/NCVOG/NCVOG.fa.tar.gz
	tar zxf NCVOG.fa.tar.gz
	rm NCVOG.fa.tar.gz
	makeblastdb -in NCVOG.fa -dbtype prot -parse_seqids
fi

IFS=,
if [[ $FORCE == 1 ]]; then
#    rm -rf VOG
echo "would remove VOG folder"
fi
if [ ! -d VOG ]; then
    mkdir -p VOG
    while read ID GENOME PROT LEN START END VOG NUM
    do
	if [ ! -z $VOG ]; then
	    blastdbcmd -db NCVOG.fa -entry $ID >> VOG/$VOG.fas
	fi
    done < NCVOG.csv
fi
make -j $CPU 
