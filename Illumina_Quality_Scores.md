#Illumina Quality Scores

# Illumina Quality Scores and FastQ files #

| Decimal | ASCII Value | Quality |
|:--------|:------------|:--------|
| 064     | `@`         | 0       |
| 065     | A           | 1       |
| 066     | B           | 2       |
| 067     | C           | 3       |
| 068     | D           | 4       |
| 069     | E           | 5       |
| 070     | F           | 6       |
| 071     | G           | 7       |
| 072     | H           | 8       |
| 073     | I           | 9       |
| 074     | J           | 10      |
| 075     | K           | 11      |
| 076     | L           | 12      |
| 077     | M           | 13      |
| 078     | N           | 14      |
| 079     | O           | 15      |
| 080     | P           | 16      |
| 081     | Q           | 17      |
| 082     | R           | 18      |
| 083     | S           | 19      |
| 084     | T           | 20      |
| 085     | U           | 21      |
| 086     | V           | 22      |
| 087     | W           | 23      |
| 088     | X           | 24      |
| 089     | Y           | 25      |
| 090     | Z           | 26      |
| 091     | `[`         | 27      |
| 092     | `\`         | 28      |
| 093     | `]`         | 29      |
| 094     | `^`         | 30      |
| 095     | `_`         | 31      |
| 096     | `````       | 32      |
| 097     | a           | 33      |
| 098     | b           | 34      |
| 099     | c           | 35      |
| 100     | d           | 36      |
| 101     | e           | 37      |
| 102     | f           | 38      |
| 103     | g           | 39      |
| 104     | h           | 40      |
| 105     | i           | 41      |
| 106     | j           | 42      |
| 107     | k           | 43      |
| 108     | l           | 44      |
| 109     | m           | 45      |
| 110     | n           | 46      |
| 111     | o           | 47      |
| 112     | p           | 48      |
| 113     | q           | 49      |
| 114     | r           | 50      |
| 115     | s           | 51      |
| 116     | t           | 52      |
| 117     | u           | 53      |
| 118     | v           | 54      |
| 119     | w           | 55      |
| 120     | x           | 56      |
| 121     | y           | 57      |
| 122     | z           | 58      |
| 123     | `{`         | 59      |
| 124     | `|`         | 60      |
| 125     | `}`         | 61      |
| 126     | `~`         | 62      |

**http://en.wikipedia.org/wiki/FASTQ_format**

Illumina 1.3+ format can encode a Phred quality score from 0 to 62 using ASCII 64 to 126 (although in raw read data Phred scores from 0 to 40 only are expected).

**From Illumina Sequencing Analysis Software User Guide For Pipeline Version 1.5:**

The quality scoring scheme has changed to the Phred scoring scheme, encoded as an ASCII character by adding 64 to the Phred value. A Phred score of a base is:
**Q<sub>phred</sub> = -10`*`log<sub>10</sub>(e)** where **e** is the estimated probability of a base being wrong.

**fastq** — This format is an adaption of the fasta format that contains quality scores.
However, the fastq format is not completely compatible with the fastq files currently
in existence, which is read by various applications (for example, BioPerl). Because a
larger dynamic range of quality scores is used, the quality scores are encoded in
ASCII as 64+score, instead of the standard 32+score. This method is used to avoid
running into non-printable characters.

The main output files in Bustard are the `_`qseq files. They have the following format:
  * 1 Machine name: (hopefully) unique identifier of the sequencer.
  * 2 Run number: (hopefully) unique number to identify the run on the sequencer.
  * 3 Lane number: positive integer (currently 1-8).
  * 4 Tile number: positive integer.
  * 5 X: x coordinate of the spot. Integer (can be negative).
  * 6 Y: y coordinate of the spot. Integer (can be negative).
  * 7 Index: positive integer. No indexing should have a value of 0.
  * 8 Read Number: 1 for single reads; 1 or 2 for paired ends.
  * 9 Sequence
  * 10 Quality: the calibrated quality string.
  * 11 Filter: Did the read pass filtering? 0 - No, 1 - Yes.