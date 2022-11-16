import cv2
import os
os.chdir("path to folder")
import random


""" docstring

    This program builds a bed file of unspliced position of de novo retrieve_transcripts
    
    input
    ------------
    - de novo transcripts
    - gtf file of transcriptome assembly
    
    output
    ------------
    - bed file
"""


def open_file(NameFile):
    F = open(NameFile, "r")
    L = F.readlines()
    return L  


def retrieve_transcripts(File):
    # This function retrieve and store transcripts from file
    L = []
    for i in File:
        if i[0] == ">":
            Name = i[1:len(i)-1]
            L.append(Name)
    return L


def make_bed_file(Liste, Info, NameFinal):
    # This function writes the bed file
    F = open(NameFinal, "w")
    for i in Info:
        ligne = i.split(",")
        TranscriptName = ligne[6]
        if TranscriptName in Liste:
            Chrom = ligne[2].split("_")[0]
            Beg = ligne[3]
            End = ligne[4]
            Direction = ligne[5]
            F.write(Chrom+"	"+Beg+"	"+End+"	"+Direction+"	"+TranscriptName+"\n")
    F.close()
    
    
def main_function():
    # Main function. Inputs have to be changed according to the line
    Fasta = open_file("AK5putativeDenovo_HighTPM.fa")
    LgoodTranscript = retrieve_transcripts(Fasta)
    Info = open_file("AK5_InfoFile")
    make_bed_file(LgoodTranscript, Info, "AK5_bed")

