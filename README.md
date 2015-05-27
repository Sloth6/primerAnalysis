# primerAnalysis

Check primer distribution in fastq files.

Usage:
	"primer_analysis.py \<primers.json\> \<myRead.fastq\>"

The primers.json must have the same format at "test_primers.json"


```
Sample output:
################################################################################
Primer analysis version 0.0.1
 	file: mysequences.fastq
	Primers: primers.json

Total reads: 10
Reads with primers: 9 (0.9)

test_primers_b
	ACAACA : 1
	ACAATA : 1
	ATAACA : 1
	ATAATA : 1

test_primers_a
	AAAAAA : 1
	GGGGGG : 1
	TTTTTT : 1

################################################################################
```