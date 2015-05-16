# wok1: https://github.com/ctb/2015-khmer-wok1-ec
# use khmer branch '2015-wok'

NULLGRAPH=../nullgraph
KHMER=../khmer
QUAKE=/Users/t/Documents/papers/2014-streaming/Quake

# rseq-mapped.fq.gz and rna.fa from 2014-streaming paper pipeline.
# ecoli-mapped.fq.gz also from 2014-streaming paper pipeline.

all: compare-sim.txt compare-rseq.txt compare-ecoli.txt \
	2015-wok-error-correction.html

clean:
	-rm simple-genome-reads.fa rseq-mapped.fq.gz.corr
	-rm ecoli-mapped.fq.gz.keep.gz

2015-wok-error-correction.html: 2015-wok-error-correction.rst
	rst2html.py 2015-wok-error-correction.rst 2015-wok-error-correction.html

simple-genome.fa:
	$(NULLGRAPH)/make-random-genome.py -l 1000 -s 1 > simple-genome.fa

simple-genome-reads.fa: simple-genome.fa
	$(NULLGRAPH)/make-reads.py -S 1 -e .01 -r 100 -C 100 simple-genome.fa --mutation-details simple-genome-reads.mut > simple-genome-reads.fa

simple-genome-reads.fa.corr: simple-genome-reads.fa
	$(KHMER)/sandbox/correct-reads.py -x 1e7 -N 4 -k 21 simple-genome-reads.fa

simple-genome-reads.sam: simple-genome-reads.fa
	bowtie2-build simple-genome.fa simple-genome > /dev/null
	samtools faidx simple-genome.fa

	bowtie2 -f -x simple-genome -U simple-genome-reads.fa -S simple-genome-reads.sam

simple-genome-sam-mismatches.pos: simple-genome-reads.sam
	./sam-scan.py simple-genome.fa simple-genome-reads.sam -o simple-genome-sam-mismatches.pos

simple-genome-reads-corr.sam: simple-genome-reads.fa.corr
	bowtie2-build simple-genome.fa simple-genome > /dev/null
	samtools faidx simple-genome.fa

	bowtie2 -f -x simple-genome -U simple-genome-reads.fa.corr -S simple-genome-reads-corr.sam

simple-genome-sam-corr-mismatches.pos: simple-genome-reads-corr.sam
	./sam-scan.py simple-genome.fa simple-genome-reads-corr.sam -o simple-genome-sam-corr-mismatches.pos

compare-sim.txt: simple-genome-sam-corr-mismatches.pos  simple-genome-sam-mismatches.pos
	./summarize-pos-file.py simple-genome-sam-mismatches.pos simple-genome-reads.fa > compare-sim.txt
	./summarize-pos-file.py simple-genome-sam-corr-mismatches.pos simple-genome-reads.fa.corr >> compare-sim.txt

############

rseq.1.bt2: rna.fa
	bowtie2-build rna.fa rseq > /dev/null
	samtools faidx rna.fa

rseq-mapped.fq.gz.corr: rseq-mapped.fq.gz
	$(KHMER)/sandbox/correct-reads.py -k 20 -Z 20 -C 3 -x 1e8 -N 4 rseq-mapped.fq.gz

rseq-corr.sam: rseq-mapped.fq.gz.corr rseq.1.bt2
	bowtie2 -p 4 -x rseq -U rseq-mapped.fq.gz.corr -S rseq-corr.sam

rseq-corr-mismatches.pos: rseq-corr.sam
	./sam-scan.py rna.fa rseq-corr.sam -o rseq-corr-mismatches.pos

rseq-mapped.sam: rseq-mapped.fq.gz rseq.1.bt2
	gunzip -c rseq-mapped.fq.gz | bowtie2 -p 4 -x rseq -U - -S rseq-mapped.sam

rseq-sam-mismatches.pos: rseq-mapped.sam
	./sam-scan.py rna.fa rseq-mapped.sam -o rseq-sam-mismatches.pos

compare-rseq.txt: rseq-corr-mismatches.pos rseq-sam-mismatches.pos
	./summarize-pos-file.py rseq-sam-mismatches.pos rseq-mapped.fq.gz > compare-rseq.txt
	./summarize-pos-file.py rseq-corr-mismatches.pos rseq-mapped.fq.gz >> compare-rseq.txt

####

ecoli-mapped.fq.gz.keep.gz: ecoli-mapped.fq.gz
	normalize-by-median.py -k 21 -x 1e8 -N 4 ecoli-mapped.fq.gz
	gzip ecoli-mapped.fq.gz.keep

ecoli-subset.100k.fq: ecoli-mapped.fq.gz
	sample-reads-randomly.py -R 2 -N 100000 ecoli-mapped.fq.gz -o ecoli-subset.100k.fq

ecoli-subset.100k.cor.fq.gz: ecoli-subset.100k.fq
	echo ecoli-subset.100k.fq > ecoli_quake_list.txt
	ln -fs /usr/local/bin/jellyfish .
	gunzip -c ecoli-mapped.fq.gz.keep.gz | $(QUAKE)/bin/count-qmers -q 33 -k 14 > ecoli_dn_counts.out
	#$(QUAKE)/bin/cov_model.py ecoli_dn_counts.out > ecoli_dn_counts.cov
	$(QUAKE)/bin/correct -f ecoli_quake_list.txt -p 4 -k 14 -q 33 -c 7.94 -z -m ecoli_dn_counts.out
	#mv ecoli-mapped.cor.fq.gz ecoli-mapped.dn.cor.fq.gz

ecoli-subset.100k.quake.fq: ecoli-subset.100k.cor.fq.gz ecoli-subset.100k.fq
	extract-original-reads-from-quake-cor.py ecoli-subset.100k.fq ecoli-subset.100k.cor.fq.gz ecoli-subset.100k.quake.fq

ecoli.dn.k21.kh: ecoli-mapped.fq.gz.keep.gz
	load-into-counting.py -k 21 -x 8e7 ecoli.dn.k21.kh ecoli-mapped.fq.gz.keep.gz

corr.k21.C5.fq: ecoli.dn.k21.kh ecoli-subset.100k.quake.fq
	$(KHMER)/sandbox/error-correct-pass2.py ecoli.dn.k21.kh ecoli-subset.100k.quake.fq -o corr.k21.C5.fq --trusted-cov 5

ecoli.1.bt2: ecoliMG1655.fa
	bowtie2-build ecoliMG1655.fa ecoli
	samtools faidx ecoliMG1655.fa

ecoli-raw.sam: ecoli.1.bt2 ecoli-subset.100k.fq
	cat ecoli-subset.100k.fq | bowtie2 -p 4 -x ecoli -U - -S ecoli-raw.sam

ecoli-raw.pos: ecoli-raw.sam
	./sam-scan.py ecoliMG1655.fa ecoli-raw.sam -o ecoli-raw.pos

ecoli-quake.sam: ecoli.1.bt2 ecoli-subset.100k.cor.fq.gz
	gunzip -c ecoli-subset.100k.cor.fq.gz | bowtie2 -p 4 -x ecoli -U - -S ecoli-quake.sam

ecoli-quake.pos: ecoli-quake.sam
	./sam-scan.py ecoliMG1655.fa ecoli-quake.sam -o ecoli-quake.pos

ecoli-corr.sam: ecoli.1.bt2 corr.k21.C5.fq
	cat corr.k21.C5.fq | bowtie2 -p 4 -x ecoli -U - -S ecoli-corr.sam

ecoli-corr.pos: ecoli-corr.sam
	./sam-scan.py ecoliMG1655.fa ecoli-corr.sam -o ecoli-corr.pos

compare-ecoli.txt: ecoli-raw.pos ecoli-quake.pos ecoli-corr.pos
	./summarize-pos-file.py ecoli-raw.pos ecoli-subset.100k.fq > compare-ecoli.txt
	./summarize-pos-file.py ecoli-quake.pos ecoli-subset.100k.cor.fq.gz >> compare-ecoli.txt
	./summarize-pos-file.py ecoli-corr.pos corr.k21.C5.fq >> compare-ecoli.txt
