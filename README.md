# Find RNAs with nontemplated repetitive tails

This workflow is designed to find nontemplated 3'-end repetitive tails in **TAIL-seq sequencing data**
```
USAGE
sbatch tail-seq.sh $1 $2 $3

COMMAND LINE VARIABLES
$1 = Sample (String. Prefix of fastq files, e.g. fastq files = ./../SRR1005384_1.fastq and ./../SRR1005384_2.fastq, Sample = SRR1005384)
$2 = Species (Possible values: human, mouse, fly, worm, worm-orsay, or arabidopsis)
$3 = Repeat (Possible values: combinations of A, T, C, G, and U)

EXAMPLE
sbatch tail-seq.sh SRR1005384 human UGUGUGUGUGUGUGUGUGUG

NOTE
Sometimes repeats are found in read names which causes the workflow to crash. Strip read names of unnecessary information by running this command: sed -i "s/ .*//" SRR1005384_2.fastq
```
