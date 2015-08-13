#Illumina QSEQ files processing

# Illumina QSEQ files processing #

According to Illumina manual **qseq** files have the following format:
  * (1) Machine name: (hopefully) unique identifier of the sequencer.
  * (2) Run number: (hopefully) unique number to identify the run on the sequencer.
  * (3) Lane number: positive integer (currently 1-8).
  * (4) Tile number: positive integer.
  * (5) X: x coordinate of the spot. Integer (can be negative).
  * (6) Y: y coordinate of the spot. Integer (can be negative).
  * (7) Index: positive integer. No indexing should have a value of 1.
  * (8) Read Number: 1 for single reads; 1 or 2 for paired ends.
  * (9) Sequence
  * (10) Quality: the calibrated quality string.
  * (11) Filter: Did the read pass filtering? 0 - No, 1 - Yes.

Tcl script **tcl\_solexa\_qseq\_parser\_2010\_01\_20.tcl** takes as an input Illumina **qseq** file and generates three output files with trimmed high quality reads:

  * **`*`.fasta** file
  * **`*`.fastq** file
  * **`*`.tab**-delimited file with structure similar to qseq format

script usage:
```
Program usage:
qseq_file_to_process, output_file, tag_length, read_length, min_length_cutoff
```

example run for data **without adaptor sequence**:
```
tclsh8.4 tcl_solexa_qseq_parser_2010_01_20.tcl  illumina_qseq_test1f_.100K.txt illumina_qseq_test1f_.100K.parsed 0 105 40
```

example run for data **with adaptor sequence**:
```
tclsh8.4 tcl_solexa_qseq_parser_2010_01_20.tcl illumina_qseq_test2A_.100K.txt illumina_qseq_test2A_.100K.parsed 7 85 47
```

see example input and output files on Download web page:
http://code.google.com/p/atgc-illumina/downloads/list

Additional trimming for remaining adaptor sequence and homogeneous tails can be done using following regular expressions:

```
-bash-3.00$ perl -p -i -e 's/AGATCGGAAGAGCGGT.*//' z-trim-test.trim
-bash-3.00$ perl -p -i -e 's/AGATCGGAAGAGCGG\n/\n/' z-trim-test.trim
-bash-3.00$ perl -p -i -e 's/AGATCGGAAGAGCG\n/\n/' z-trim-test.trim
-bash-3.00$ perl -p -i -e 's/AGATCGGAAGAGC\n/\n/' z-trim-test.trim
-bash-3.00$ perl -p -i -e 's/AGATCGGAAGAG\n/\n/' z-trim-test.trim
-bash-3.00$ perl -p -i -e 's/AGATCGGAAGA\n/\n/' z-trim-test.trim
-bash-3.00$ perl -p -i -e 's/AGATCGGAAG\n/\n/' z-trim-test.trim
-bash-3.00$ perl -p -i -e 's/AGATCGGAA\n/\n/' z-trim-test.trim
-bash-3.00$ perl -p -i -e 's/AGATCGGA\n/\n/' z-trim-test.trim

-bash-3.00$ perl -p -i -e 's/^A{8,}//' z-trim-test.trim
-bash-3.00$ perl -p -i -e 's/^T{8,}//' z-trim-test.trim
-bash-3.00$ perl -p -i -e 's/^C{8,}//' z-trim-test.trim
-bash-3.00$ perl -p -i -e 's/^G{8,}//' z-trim-test.trim

-bash-3.00$ perl -p -i -e 's/A{8,}$//' z-trim-test.trim
-bash-3.00$ perl -p -i -e 's/T{8,}$//' z-trim-test.trim
-bash-3.00$ perl -p -i -e 's/G{8,}$//' z-trim-test.trim
-bash-3.00$ perl -p -i -e 's/C{8,}$//' z-trim-test.trim
```