import cv2
import os
os.chdir("path to folder")
import random


""" docstring

    This program provides statistics about how many positions are filled by de novo transcripts
    
    input
    ------------
    - output from categories attribution (Pg7)
    
    output
    ------------
    - stats on positions
"""


def open_file(NameFile):
    F = open(NameFile, "r")
    L = F.readlines()
    return L  


def make_dic_position(FilePositions):
    # This function retrieve the genomic position of transcripts and store them in a dictionary
    Dico = {}
    for i in FilePositions:
        ligne = i.split("\n")[0]
        ligne = ligne.split(",")
        Dico[ligne[0]] = ligne[1]
    return Dico


def count_position(Dico):
    # This function builds the final statistics
    UndirectionalExonLonger = 0
    UndirectionalExonInside = 0
    Intergenic = 0
    ReverseGenic = 0
    Pseudogene = 0
    NcRNA = 0
    Intronic = 0
    Weird = 0
    ExonLonger = 0
    ExonInside = 0
    for i in Dico.keys():
        PosORF = Dico[i]
        if PosORF == "UndirectionalExonLonger":
            UndirectionalExonLonger+=1
        elif PosORF == "UndirectionalExonInside":
            UndirectionalExonInside+=1
        elif PosORF == "Intergenic":
            Intergenic+=1
        elif PosORF == "ReverseGenic":
            ReverseGenic+=1
        elif PosORF == "Pseudogene":
            Pseudogene+=1
        elif PosORF == "NcRNA":
            NcRNA+=1
        elif PosORF == "ExonLonger":
            ExonLonger+=1
        elif PosORF == "ExonInside":
            ExonInside+=1
        elif PosORF == "Intronic":
            Intronic+=1
        else:
            print PosORF
            Weird+=1
    print ("UndirectionalExonLonger : "+str(UndirectionalExonLonger))
    print ("UndirectionalExonInside : "+str(UndirectionalExonInside))
    print ("ExonInside : "+str(ExonInside))
    print ("ExonLonger : "+str(ExonLonger))
    print ("Intergenic : "+str(Intergenic))
    print ("ReverseGenic : "+str(ReverseGenic))
    print ("Pseudogene : "+str(Pseudogene))
    print ("NcRNA : "+str(NcRNA))
    print ("Intronic : "+str(Intronic))
    print ("Weird : "+str(Weird))
    

def main_function():
    # Has to be modified acording to the line
    FilePositions = open_file("AK5_PositionsInGenome")
    Dico = make_dic_position(FilePositions)
    count_position(Dico)
