#!/bin/bash
#$ -S /bin/bash # tell the SGE to use bash
#$ -cwd         #run in current working dir
#$ -pe smp 4    # how many CPUs?
#$ -l h_vmem=3G # how much RAM?
#$ -w e         # reject jobs with error
#$ -V           # export all environment variables
#$ -N GetIntrons # job name



perl extract_seq_from_gff3.pl -d TotaleGenome.fa - Drosophila_melanogaster.BDGP6.32.54.chr.gff3 > output_intron.fa
