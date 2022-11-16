import cv2
import os
os.chdir("path to folder")
import random


""" docstring

    This program builds fasta file according to the genomic position of the transcript
    
    input
    ------------
    - FASTA file of transcripts
    
    output
    ------------
    - several FASTA files of transcripts at the same genomic position
    
"""


def open_file(NameFile):
    F=open(NameFile, "r")
    L=F.readlines()
    return L  


def build_dic_position(ListePositions, Header):
    # This function builds a dictionary containing the transcripts of each positions in genome
    D = {}
    for i in ListePositions:
        ligne = i.split("\n")[0]
        ligne = ligne.split(",")
        Name = Header+"_"+ligne[0]
        Cat = ligne[1]
        D [Name] = Cat
    return D


def prepare_fasta_files(Fasta, Dico, Header):
    # This function builds the fasta files 
    name_f1 = Header+"_"+"NcRNA"
    name_f2 = Header+"_"+"ExonLonger"
    name_f3 = Header+"_"+"Intergenic"
    name_f4 = Header+"_"+"ReverseGenic"
    name_f5 = Header+"_"+"Pseudogene"
    name_f6 = Header+"_"+"Intronic"
    F1 = open(name_f1, "w")
    F2 = open(name_f2, "w")
    F3 = open(name_f3, "w")
    F4 = open(name_f4, "w")
    F5 = open(name_f5, "w")
    F6 = open(name_f6, "w")
    for i in Fasta:
        if i[0] == ">":
            CompleteName = i
            Name = i[1:len(i)-1]
            Cat = Dico[Name]
        else:
            Seq = i
            if Cat == "NcRNA":
                F1.write(CompleteName)
                F1.write(Seq)
            elif Cat == "ExonLonger":
                F2.write(CompleteName)
                F2.write(Seq)
            elif Cat == "Intergenic":
                F3.write(CompleteName)
                F3.write(Seq)
            elif Cat == "ReverseGenic":
                F4.write(CompleteName)
                F4.write(Seq)
            elif Cat == "Pseudogene":
                F5.write(CompleteName)
                F5.write(Seq)
            elif Cat == "Intronic":
                F6.write(CompleteName)
                F6.write(Seq)
    F1.close()
    F2.close()
    F3.close()
    F4.close()
    F5.close()
    F6.close()
            

def main_function():
    # Main function. The input have to be modified according to the line of interest
    ListePositions = open_file("AK5_PositionsInGenome")
    Dico = build_dic_position(ListePositions, "AK5")
    ListeFasta = open_file("AK5_CleanedDeNovo_Renamed.fa")
    prepare_fasta_files(ListeFasta, Dico, "AK5")
