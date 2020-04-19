

seqs = $(wildcard *.fas)
cons = $(seqs:.fas=.cons)

all: VOG.consensus.fas

.PHONY: clean
clean:
	rm -f *.hmm *.cons *.fasaln

VOG.consensus.fas: $(cons)
	echo $(cons)
	cat $(cons) > ../VOG.consensus.fas
	
%.fasaln : %.fas
	muscle -in $< -out $@

%.hmm : %.fasaln
	hmmbuild $@ $< 

%.cons : %.hmm
	hmmemit -C -o $@ $<
