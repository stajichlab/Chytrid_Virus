

seqs = $(wildcard Alignments/Post-trimming/*.faa)
cons = $(seqs:.faa=.cons)
hmm  = $(seqs:.faa=.hmm)

all: NCLDV.consensus.fasta NCLDV.hmm NCLDV.hmmb

.PHONY: clean
clean:
	rm -f *.hmm *.cons

%.hmmb: %.hmm
	hmmconvert -b $< > $@
	hmmpress $@

NCLDV.hmm: $(hmm)
	cat $(hmm) > $@
	hmmpress $@

NCLDV.consensus.fasta: $(cons)
	cat $(cons) > $@
	
%.hmm : %.faa
	hmmbuild $@ $< 

%.cons : %.hmm
	hmmemit -C -o $@ $<
