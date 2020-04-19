#!/usr/bin/bash
#SBATCH -p short

if [ ! -d spades_noreadfilter ]; then
	# replace with a osf.io download in future?
	gdrive download -r 1Bv1MkQzki9GJLpE8iziCSaUuhD5hk1_s
	pigz -d spades_noreadfilter/*.gz
fi
pushd db/NCVOG

# next make the VOG files
# or replace with
#bash make_NCVOG_alns.sh
sbatch make_NCVOG_alns.sh

