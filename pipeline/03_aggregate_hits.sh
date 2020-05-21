#!/usr/bin/bash
#SBATCH -p short --mem 8gb --out logs/aggregate.log -N 1 -n 6
module unload perl
module load parallel
cat | parallel -j 6  << EOF
./scripts/hmmsearch_table.py --prediction prodigal --library NCLDV
./scripts/hmmsearch_table.py --prediction prodigal --library NCVOG 
./scripts/hmmsearch_table.py --prediction prodigal --library NCLDV --nofilter 
./scripts/hmmsearch_table.py --prediction prodigal --library NCVOG --nofilter
./scripts/hmmsearch_table.py --prediction funannotate --library NCLDV
./scripts/hmmsearch_table.py --prediction funannotate --library NCVOG
./scripts/NCVOGnegctlsearch_table.py
EOF

#./scripts/NCVOGsearch_table.py
