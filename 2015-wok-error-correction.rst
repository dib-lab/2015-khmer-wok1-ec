

One of the features in khmer that we're pretty excited about is
the read-to-graph aligner, which gives us a way to align sequences
to a De Bruijn graph; our nickname for it is "graphalign."

Briefly, graphalign uses a pair-HMM to align a sequence to a k-mer
graph (aka De Bruijn graph) allowing both mismatches and indels, and
taking into account coverage using a binary model (trusted and
untrusted k-mers).  The core code was written by Jordan Fish when he
was a graduate student in the lab, and was then refactored by Michael
Crusoe.

Graphalign actually lets us do lots of things, including align both
short and long sequences to DBG graphs, error correct, and call
variants.  We've got a simple API built into khmer, and we're working
to build it out.

The core API is based around the concept of a ReadAligner object::

    aligner = khmer.ReadAligner(graph, trusted_cov, bits_theta)

where 'trusted_cov' defines what the trusted k-mer coverage is, and
bits_theta adjusts a scoring parameter used to extend alignments.

This object can be used to align short sequences to the graph::

     score, graph_alignment, read_alignment, truncated = \
         aligner.align(read)

Here, 'graph_alignment' and 'read_alignment' are strings; if
'truncated' is false, then they are of the same length, and constitute
a full gapped alignment of the DNA sequence in 'read' to the graph.

The approach used by 'align' is to seed an alignment at the first trusted
k-mer, and then extend the alignment along the graph in both directions.
Thus, it's effectively a local aligner.

Error correction
~~~~~~~~~~~~~~~~

Our initial motivation for graphalign was to use it to do error
correction, with specific application to short-read sequences.  There
was (and to some extent still is) a dearth of error correction
approaches that can be used for metagenome and transcriptome data
sets, and since that kind of data is what our lab works on, we needed
an error correction approach for those data.  We also wanted something
a bit more programmable than the existing error correctors, which were
primarily command-line tools; we've found a lot of value in building
libraries, and wanted to use that approach here, too.

The basic idea is this: we build a graph from our short-read data,
and then go back through and align each short read to the graph.  A
successful alignment is then the corrected read.  The basic code looks
like this::

    graph = build_graph(dataset)

    aligner = khmer.ReadAligner(graph, trusted_cov, bits_theta)

    for read in dataset:
        score, ga, ra, is_trunc = aligner.align(read)
        corrected_read = ga

In conjunction with `our work on semi-streaming algorithms
<https://peerj.com/preprints/890/>`__, we can directly convert this
into a semi-streaming algorithm that works on genomes, metagenomes,
and transcriptomes.  This is implemented in the `correct-errors script
<https://github.com/dib-lab/khmer/blob/2015-wok/sandbox/correct-reads.py>`__.

Some results
~~~~~~~~~~~~

If we try this out on a simulated data set (random genome, randomly
chosen reads - see target ``compare-sim.txt`` in `the Makefile
<https://github.com/ctb/2015-khmer-wok1-ec/blob/master/Makefile>`__),
it takes the simulated data from an error rate of around 1% to about
0.1%; see `compare-sim.txt
<https://github.com/ctb/2015-khmer-wok1-ec/blob/master/compare-sim.txt>`__.

Appying this to the same ~750k read subset of mRNAseq that we tackled
in the semi-streaming paper (the data itself is from the `Trinity
paper, Grabherr et al, 2011
<http://www.ncbi.nlm.nih.gov/pubmed/21572440>`__), we take the data
from an error rate of about 2.1% to 1.2% (see target
``rseq-compare.txt`` in `Makefile
<https://github.com/ctb/2015-khmer-wok1-ec/blob/master/Makefile>`__).
There are several reasons why this doesn't correct many errors -
first, error correction depends on high coverage, and much of this
RNAseq data set is low coverage; second, this data set has a lot of
errors; and third, RNAseq may have a broader k-mer abundance
distribution than genomic sequencing.

One important side note: we use exactly the same script for error
correcting RNAseq data as use do for genomic data.

How good is the error correction?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

tl; dr? It's pretty good but not better than current methods.  When we
compare to Quake results on an E. coli data set (target
`compare-ecoli.txt
<https://github.com/ctb/2015-khmer-wok1-ec/blob/master/compare-ecoli.txt>`__
in `the Makefile
<https://github.com/ctb/2015-khmer-wok1-ec/blob/master/Makefile>`__),
we see:

============  ==========
Data set      Error rate
============  ==========
Uncorrected   1.576%
Quake         0.010%
khmer         0.015%
============  ==========

(Note that here we do a fair comparison by looking only at errors on
sequences that Quake doesn't discard.)

Concluding thoughts
~~~~~~~~~~~~~~~~~~~

What attracts us to this approach is that it's really *simple*.
