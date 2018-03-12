Introduction
============

AC-DIAMOND attempts to speed up `DIAMOND <http://github.com/bbuchfink/diamond>`_ via better SIMD parallelization and compressed indexing. Experimental results show that AC-DIAMOND was about 6~7 times faster than DIAMOND on aligning DNA reads or contigs while retaining the essentially the similar sensitivity. AC-DIAMOND was developped based on DIAMOND v0.7.9.

Geting Started
==============
::

    Download the source code and install boost: ./install-boost
    Compile the source code: make
    ./ac-diamond makedb --in nr.faa -d nr
    ./ac-diamond align -d nr -q reads.fna -a matches -t <temporary directory>
    ./ac-diamond view -a matches.daa -o matches.m8

Commands
========
Commands are issued as the first parameter on the command line and set the task to be run by the program.

======= ===========
Command Description
======= ===========
makedb  Create AC-DIAMOND formatted reference database from a FASTA input file.
align  Align translated DNA query sequences against a protein reference database.
view    Generate formatted output from DAA files.
======= ===========

Makedb options
==============
============ ===== ======= ===========
Option       Short Default Description
============ ===== ======= ===========
--threads    -p    max     Number of CPU threads.
--in                       Path to protein reference database file in FASTA format (may be gzip compressed).
--db         -d            Path to AC-DIAMOND database file.
--block-size -b    4       Block size in billions of sequence letters to be processed at a time.
--sensitive                Build index for searching alignments in the sensitive mode.
============ ===== ======= ===========

General options
====================
=================== ===== ======= ===========
Option              Short Default Description
=================== ===== ======= ===========
--threads           -p    max     Number of CPU threads.
--db                -d            Path to AC-DIAMOND database file (not including the file extension).
--query             -q            Path to query input file in FASTA or FASTQ format (may be gzip compressed).
--query-block-size  -z    6       query sequence block size in billions of letters
--daa               -a            Path to output file in DAA format (extension .daa will be appended).
=================== ===== ======= ===========

Sensitivity & speed options
===========================
=========== ===== ======= ===========
Option      Short Default Description
=========== ===== ======= ===========
--sensitive               Trigger the sensitive alignment mode with a 10x8 seed shape configuration.
--band            auto    Dynamic programming band for seed extension. This corresponds to the maximum length of gaps that can be found in alignments.
=========== ===== ======= ===========

Scoring & Reporting Options
===========================
================= ===== ======== ===========
Option            Short Default  Description
================= ===== ======== ===========
--gapopen               11       Gap open penalty.
--gapextend             1        Gap extension penalty.
--matrix                BLOSUM62 Scoring matrix.
--seg                            Enable SEG masking of low complexity segments in the query (yes/no). The default is no for blastp and yes for blastx.
--max-target-seqs -k    25       The maximum number of target sequences per query to keep alignments for.
--top                            Keep alignments within the given percentage range of the top alignment score for a query (overrides –max-target-seqs option).
--evalue          -e    0.001    Maximum expected value to keep an alignment.
--min-score                      Minimum bit score to keep an alignment. Setting this option will override the --evalue parameter.
================= ===== ======== ===========

Memory & performance options
============================
============== ===== ======== ===========
Option         Short Default  Description
============== ===== ======== ===========
--tmpdir       -t    /dev/shm Directory to be used for temporary storage.
============== ===== ======== ===========

View options
============
========== ===== ======== ===========
Option     Short Default  Description
========== ===== ======== ===========
--daa      -a             Path to input file in DAA format.
--out      -o             Path to output file.
--outfmt   -f             Format of output file. (tab = BLAST tabular format; sam = SAM format)
--compress       0        Compression for output file (0=none, 1=gzip).
========== ===== ======== ===========
