#import cv2
import os
from Bio.Seq import Seq
os.chdir("path to folder")
import random


""" docstring

    This program extract TPM values of de novo transcripts, and build a list of de novo transcript with TPM value higher than 0.5
    
    input
    ------------
    - gtf file
    - de novo transcripts
    
    output
    ------------
    - final de novo transcripts
"""


def open_file(NameFile):
    F = open(NameFile, "r")
    L = F.readlines()
    return L  


def build_dic_TPM(File):
    # This function builds a dictionary of transcripts associated to their tpm value
    Dico = {}
    for i in File:
        ligne = i.split("\n")[0]
        ligne = ligne.split(",")
        TotalName = ligne[6]
        TPM = ligne[13]
        Dico[TotalName] = TPM
    return Dico


def build_final_file(DicoTPM, FastaFile, FinalFileName):
    # This function builds the final output file
    NbTotalSequences=0
    DicoSeqs = {}
    for i in FastaFile:
        if i[0] == ">":
            LongName = i
            Name  =i[1:len(i)-1]
            TPM = float(DicoTPM[Name])
            NbTotalSequences+=1
        else:
            Seq = i
            if TPM>=0.5:
                DicoSeqs[LongName] = Seq
    print ("Total Number of putative de novo transcripts : " + str(NbTotalSequences))
    print ("Total Number of putative de novo transcripts with RPM>0,5 : " + str(len(DicoSeqs)))
    F = open(FinalFileName, "w")
    for i in DicoSeqs:
        F.write(i)
        F.write(DicoSeqs[i])
    F.close()


def main_function():
    # Main function. Has to be modified according to the name of the line
    FastaFile = open_file("AK5_FASTA_DeNovo_Trancript_Forward_NonDipteran")
    InfoFile = open_file("AK5_InfoFile")
    DicoTPM = build_dic_TPM(InfoFile)
    build_final_file(DicoTPM, FastaFile, "AK5putativeDenovo_HighTPM.fa")
