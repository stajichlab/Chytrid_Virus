

seqs = $(wildcard VOG/*.fas)
cons = $(seqs:.fas=.cons)
hmm  = $(seqs:.fas=.hmm)

all: VOG.consensus.fas VOG.hmm

.PHONY: clean
clean:
	rm -f *.hmm *.cons *.fasaln

VOG.hmm: $(hmm)
	cat $(hmm) > $@

VOG.consensus.fas: $(cons)
	cat $(cons) > $@
	
%.fasaln : %.fas
	muscle -in $< -out $@

%.hmm : %.fasaln
	hmmbuild $@ $< 

%.cons : %.hmm
	hmmemit -C -o $@ $<
