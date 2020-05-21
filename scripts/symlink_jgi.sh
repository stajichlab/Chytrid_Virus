#!/usr/bin/bash

pushd proteins
for file in $(ls *.aa.fasta | grep -v LCG); do 
	f=$(basename $file .aa.fasta); 
	if [ -s ../gff3/$f.JGI.genes.gff3 ]; then
		continue
	fi
	if [ -f /bigdata/stajichlab/shared/projects/1KFG/2019_dataset/data/GFF3/$f.gff3.gz ]; then 
		echo $file;  
		zgrep -P "\tgene\t" /bigdata/stajichlab/shared/projects/1KFG/2019_dataset/data/GFF3/$f.gff3.gz > ../gff3/$f.JGI.genes.gff3
	else
		echo "cannot find $f.gff3.gz in JGI folder"	
	fi
done
popd
