import os
os.chdir("path to folder")
import random


""" docstring

    This program extract all transcript in transcriptome assembly of a line and rename them
    
    input
    ------------
    - transcriptome assembly (FASTA)
    
    output
    ------------
    -new name of transcripts in transcriptome assembly
"""


def open_file(NameFile):
    F=open(NameFile, "r")
    L=F.readlines()
    return L  


def create_seqfile_new_names(File):
    renamed_transcrips_pop = open("AK5RenamedTranscripts", "w")
    for i in File:
        if i[0] == ">":
            ligne = i.split(" ")
            NewName = ligne[0]+"_"+ligne[1]
            renamed_transcrips_pop.write(NewName)
        else:
            renamed_transcrips_pop.write(i)
    renamed_transcrips_pop.close()


def main_function():
    File = open_file("AK5assemblyTranscripts.fa")
    create_seqfile_new_names(File)
