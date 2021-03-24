# mFinder

This set of scripts was developed to identify and to generate intrahost/co-infection consensus genomes of SARS-CoV-2. The test sequences were generate by reference genome assembly with bowtie2 usign the SRR13418742 data from Gupta et al. 2021 [REF](https://www.microbiologyresearch.org/content/journal/jgv/10.1099/jgv.0.001562).

## Preparing the Enviroment

This script need 2 pre-installed tools:

- [bam-readcount](https://github.com/genome/bam-readcount)
- [MAFFT](https://mafft.cbrc.jp/alignment/software/)

With the bam-readcount and MAFFT intalled, you can use conda to install all dependencies:

> conda env create -f mFinder.yml

> conda activate mFinder

## Usage
The main.sh compile all steps to perform the analysis, and 3 files should be parsed:
- reference = reference wuhan genome (NC_045512.2)
- consensus = the consensus genome that will be screened to intrahost variants
- sorted_bam = sorted bam file that was used to generate the consensus sequences


- Running 
> bash main.sh reference consensus sorted_bam

- Test files:
> bash main.sh test_files/NC_045512v2.fa test_files/SRR13418742.consensus.fa test_files/SRR13418742.sorted.bam


## Outputs
- *.tsv.fmt.minors.tsv.fmt - tab separed file with the nucleotide diversity by intrahost variant postition, the nucleotide and depth of minor and major variants
- *.algn.major.fa - fasta alignment with the major variant consensus sequence
- *.algn.minor.fa - fasta alignment with the minor variant consensus sequence

## Some limitations
- The name of consensus and sorted_bam files should contain the sequence name header of consenus sequences;
- This set of scripts didn't recovery intrahost variants related do indel regions, these region ones should be investigated manually.

## Disclaimer

- This script will continue to be developed to avoid erros related do sequence name, to remove the reference genome of output, to work with multiple sequences and to work with indel regions.