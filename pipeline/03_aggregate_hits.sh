#!/usr/bin/bash
#SBATCH -p short

./scripts/NCLDVhmmsearch_table.py
./scripts/NCVOGhmmsearch_table.py
./scripts/NCVOGnegctlsearch_table.py
./scripts/NCVOGsearch_table.py
