#!/bin/bash
#$ -S /bin/bash # tell the SGE to use bash
#$ -cwd         #run in current working dir
#$ -pe smp 16    # how many CPUs?
#$ -l h_vmem=3G # how much RAM?
#$ -w e         # reject jobs with error
#$ -V           # export all environment variables
#$ -N AK5_Prots     # job name





gffread -g AK5finalGenome.masked.fa -y Proteins.fa AK5deNovoGenomeAnnotation.gff


