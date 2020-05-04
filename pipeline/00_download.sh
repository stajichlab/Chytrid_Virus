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
