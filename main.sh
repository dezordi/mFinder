#!/bin/bash
#Written by: Filipe Dezordi (https://dezordi.github.io/)
#At FioCruz/IAM - 11 Nov. 2020


##USAGE:

#Arguments
REFERENCE=$1 #reference genome
FASTA=$2 #consensus SARS-CoV-2 genome
BAM=$3 #sorted bam file used to generate the consensus SARS-CoV-2 genome

bam-readcount -d 50000 -b 30 -q 30 -w 0 -f $REFERENCE $BAM > $BAM.bamreadcount.tsv
python minor_finder.py -in $BAM.bamreadcount.tsv
python major_minor.py -in $BAM.bamreadcount.tsv.fmt.minors.tsv
mafft --thread 1 --keeplength --add $FASTA $REFERENCE > $FASTA.ref.algn
python put_minor.py -in $FASTA.ref.algn -mv $BAM.bamreadcount.tsv.fmt.minors.tsv.fmt