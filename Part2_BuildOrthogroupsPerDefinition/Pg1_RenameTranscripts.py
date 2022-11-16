import cv2
import os
os.chdir("path to folder")
import random


""" docstring

    This program modify the name of de novo transcripts
    
    input
    ------------
    - FASTA file of de novo transcripts
    
    output
    ------------
    - FASTA file of de novo transcripts renamed
    
"""


def open_file(NameFile):
    F = open(NameFile, "r")
    L = F.readlines()
    return L  


def modify_file(File, Header, FinalName):
    F = open(FinalName, "w")
    for i in File:
        if i[0] == ">":
            NewName = i[0]+Header+"_"+i[1:]
            F.write(NewName)
        else:
            F.write(i)
    F.close()


def main_function():
    AK5 = open_file("AK5_CleanedDeNovo.fa")
    modify_file(AK5, "AK5", "AK5_CleanedDeNovo_Renamed.fa")
    DK5 = open_file("DK5_CleanedDeNovo.fa")
    modify_file(DK5, "DK5", "DK5_CleanedDeNovo_Renamed.fa")
    GI5 = open_file("GI5_CleanedDeNovo.fa")
    modify_file(GI5, "GI5", "GI5_CleanedDeNovo_Renamed.fa")
    SW5 = open_file("SW5_CleanedDeNovo.fa")
    modify_file(SW5, "SW5", "SW5_CleanedDeNovo_Renamed.fa")
    UM = open_file("UM_CleanedDeNovo.fa")
    modify_file(UM, "UM", "UM_CleanedDeNovo_Renamed.fa")
    YE = open_file("YE_CleanedDeNovo.fa")
    modify_file(YE, "YE", "YE_CleanedDeNovo_Renamed.fa")
    Zamb = open_file("Zamb_CleanedDeNovo.fa")
    modify_file(Zamb, "Zamb", "Zamb_CleanedDeNovo_Renamed.fa")
