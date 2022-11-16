import cv2
import os
from Bio.Seq import Seq
os.chdir("path to folder")
import random


""" docstring

    This program build a FASTA file with only transcripts that have no BLAST hits
    
    input
    ------------
    - output file from pg3 with name of de novo transcripts
    - renamed transcriptome assembly
    
    output
    ------------
    - FASTA file of de novo transcripts
"""

def open_file(NameFile):
    F = open(NameFile, "r")
    L = F.readlines()
    return L  


def make_final_file(DeNovo, AllTransc, NameFile):
    Liste = []
    Dico  = {}
    for i in DeNovo:
        ligne = i.split("\n")[0]
        Liste.append(ligne)
    for i in AllTransc:
        if i[0] == ">":
            Name = i.split("\n")[0][1:]
        else:
            Transcript = i.split("\n")[0]
            Dico[Name] = Transcript
    F = open(NameFile, "w")
    for i in Liste:
        if i in Dico.keys():
            F.write(">"+i+"\n")
            F.write(Dico[i]+"\n")
    F.close()


main_function():
    # has to be modified according to th line
    BigFile = open_file("path to folder/AK5_FASTA_DeNovoTranscForward")
    DeNovo = open_file("AK5_DeNovo_Trancript_Forward_NonDipteran")
    make_final_file(DeNovo, BigFile, "AK5_FASTA_DeNovo_Trancript_Forward_NonDipteran")

