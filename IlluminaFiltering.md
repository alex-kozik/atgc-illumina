#Introduction

# Introduction #

BaseCalling

SCRIPT: tcl\_solexa\_sig2\_base\_caller\_023\_200-100-60-1KQ.tcl - Tcl script for Solexa basecalling

```

bash-2.03$
bash-2.03$ tclsh8.3 tcl_solexa_sig2_base_caller_023_200-100-60-1KQ.tcl
Program usage:
Table_to_Process, output_file_good, output_file_bad, output_file_base_call, id_prefix
bash-2.03$

bash-2.03$ tclsh8.4 tcl_solexa_sig2_base_caller_023_200-100-60-1KQ.tcl solexa_sig2_example.txt solexa_sig2_example.basecall.good solexa_sig2_example.basecall.bad solexa_sig2_example.basecall.debug.tab SEQ_TEST3_

```

INPUT: solexa\_sig2\_example.txt - example 'sig2' input file

OUTPUT:
  * solexa\_sig2\_example.basecall.bad
  * solexa\_sig2\_example.basecall.debug.tab
  * solexa\_sig2\_example.basecall.good
  * solexa\_sig2\_example.basecall.good.qual
  * solexa\_sig2\_example.basecall.good.trim.fasta
  * solexa\_sig2\_example.basecall.good.trim.qual

Note: quality scores are in following scale:
A - excellent
B - very good
C - good
D - acceptable
F - not reliable

see the latest version of the script at http://code.google.com/p/atgc-illumina/source/browse/#svn/trunk