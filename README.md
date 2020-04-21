This is a project for finding Viral DNA genomes from assembled Chytrid projects.

Collaboration with Jason Stajich Lab (http://lab.stajich.org) and Tim James (http://umich.edu/~mycology/people.html)

Download
==
To download the JGI dataset you need to copy the `scripts/init_jgi_download.sh.template` to `scripts/init_jgi_download.sh` and put in your own username and password for download access.

Then run 

```
bash scripts/init_jgi_download.sh
```

This depends on python3 which runs with some XML parsing but it should be core python.

Running
==
0. Setup the NCVOG db - go into db/NCVOG and run `make_NCVOG_alns.sh` - requires HMMER3, NCBI-BLAST, muscle, and parallel 
1. Translated searches - requires FASTA36 installed
2. HMM searches - requires HMMER3 installed
