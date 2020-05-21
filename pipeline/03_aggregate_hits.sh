#!/usr/bin/bash
#SBATCH -p short

./scripts/hmmsearch_table.py --prediction prodigal --library NCLDV
./scripts/hmmsearch_table.py --prediction prodigal --library NCVOG
./scripts/hmmsearch_table.py --prediction prodigal --library NCLDV --nofilter
./scripts/hmmsearch_table.py --prediction prodigal --library NCVOG --nofilter

./scripts/hmmsearch_table.py --prediction fuannotate --library NCLDV
./scripts/hmmsearch_table.py --prediction funannotate --library NCVOG

./scripts/NCVOGnegctlsearch_table.py
#./scripts/NCVOGsearch_table.py
