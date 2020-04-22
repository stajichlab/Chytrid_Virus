

seqs = $(wildcard NCVOG/*.fas)
cons = $(seqs:.fas=.cons)
hmm  = $(seqs:.fas=.hmm)

all: NCVOG.consensus.fasta NCVOG.hmm NCVOG.hmmb

.PHONY: clean
clean:
	rm -f *.hmm *.cons *.fasaln

%.hmmb: %.hmm
	hmmconvert -b $< > $@
	hmmpress $@

NCVOG.hmm: $(hmm)
	cat $(hmm) > $@
	hmmpress $@

NCVOG.consensus.fasta: $(cons)
	cat $(cons) > $@
	
%.fasaln : %.fas
	muscle -in $< -out $@

%.hmm : %.fasaln
	hmmbuild $@ $< 

%.cons : %.hmm
	hmmemit -C -o $@ $<
