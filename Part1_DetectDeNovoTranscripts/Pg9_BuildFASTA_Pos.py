import cv2
import os
os.chdir("path to folder")
import random


""" docstring

    This program build separated FASTA files of de novo transcripts based on their genomic AK5_PositionsInGenome
    
    input
    ------------
    - position in genomes
    - FASTA file of de novo transcripts
    
    output
    ------------
    - Separated FASTA files
"""

def open_file(NameFile):
    F=open(NameFile, "r")
    L=F.readlines()
    return L  


def build_list_cat(Names):
    L = []
    for i in Names:
        ligne = i.split(",")
        #print ligne
        Name = ligne[0]
        Cat = ligne[1][0:len(ligne[1])-1]
        if Cat != "UndirectionalExonInside" and Cat != "ExonInside":
            L.append(Name)
    #print len(L)
    return L


def build_final_file(Lgoods, FASTA, NameFile):
    F = open(NameFile, "w")
    for i in FASTA:
        if i[0] == ">":
            Name = i[1:len(i)-1]
            LongName = i
        else:
            Seq = i
            if Name in Lgoods:
                F.write(LongName)
                F.write(Seq)
    F.close()


def main_function():
    Names = open_file("AK5_PositionsInGenome")
    FASTA = open_file("AK5putativeDenovo_HighTPM.fa")
    Lgoods = build_list_cat(Names)
    build_final_file(Lgoods, FASTA, "AK5_CleanedDeNovo.fa")
