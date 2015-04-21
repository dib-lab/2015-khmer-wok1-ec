NULLGRAPH=../nullgraph
KHMER=../khmer

all: simple-genome-reads.fa simple-genome-reads.fa.corr \
	simple-genome-sam-corr-mismatches.pos \
	simple-genome-sam-mismatches.pos compare-sim.txt

clean:
	-rm simple-genome-reads.fa

simple-genome.fa:
	$(NULLGRAPH)/make-random-genome.py -l 1000 -s 1 > simple-genome.fa

simple-genome-reads.fa: simple-genome.fa
	$(NULLGRAPH)/make-reads.py -S 1 -e .01 -r 100 -C 100 simple-genome.fa --mutation-details simple-genome-reads.mut > simple-genome-reads.fa

simple-genome-reads.fa.corr: simple-genome-reads.fa
	$(KHMER)/sandbox/correct-errors.py -x 1e7 -N 4 -k 20 simple-genome-reads.fa

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
