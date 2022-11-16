import cv2
import os
os.chdir("path to folder")
import random


""" This program merge together the info file of each de novo transcripts of each ligne

    input
    ------------
    - info files of each ligne
    
    output
    ------------
    - All de novo transcripts info files
"""


def open_file(NameFile):
    F=open(NameFile, "r")
    L=F.readlines()
    return L  


def dic_position(File, Header):
    Dico = {}
    for i in File:
        ligne = i.split(",")
        Name = Header+"_"+ligne[0]
        Dico[Name] = ligne[1] 
    return Dico


def update_info(File, Dico):
    Liste = []
    for i in File[1:]:
        ligne = i.split("\n")[0]
        lolo = ligne.split(",")
        if lolo[6] in Dico.keys():
            NewLigne = ligne+","+Dico[lolo[6]]
            Liste.append(NewLigne)
    return Liste


def build_final_file(L1, L2, L3, L4, L5, L6, L7, Name):
    F = open(Name, "w")
    F.write("GeneName,NbTranscripts,Chrom,PosBeg,PosEnd,Direction,TotalTranscName,TranscriptName,UnsplicedSize,SplicedSize,NbExon,Cov,FPKM,TPM,Category"+"\n")
    for i in L1:
        F.write(i)
    for i in L2:
        F.write(i)
    for i in L3:
        F.write(i)
    for i in L4:
        F.write(i)
    for i in L5:
        F.write(i)
    for i in L6:
        F.write(i)
    for i in L7:
        F.write(i)
    F.close()
        

def main_function():
    Position_AK5 = open_file("AK5_PositionsInGenome")
    Dico_AK5 = dic_position(Position_AK5, "AK5")
    Info_AK5 = open_file("AK5_New_InfoFile")
    ListeAK5 = update_info(Info_AK5, Dico_AK5)
    Position_DK5 = open_file("DK5_PositionsInGenome")
    Dico_DK5 = dic_position(Position_DK5, "DK5")
    Info_DK5 = open_file("DK5_New_InfoFile")
    ListeDK5 = update_info(Info_DK5, Dico_DK5)
    Position_GI5 = open_file("GI5_PositionsInGenome")
    Dico_GI5 = dic_position(Position_GI5, "GI5")
    Info_GI5 = open_file("GI5_New_InfoFile")
    ListeGI5 = update_info(Info_GI5, Dico_GI5)
    Position_SW5 = open_file("SW5_PositionsInGenome")
    Dico_SW5 = dic_position(Position_SW5, "SW5")
    Info_SW5 = open_file("SW5_New_InfoFile")
    ListeSW5 = update_info(Info_SW5, Dico_SW5)
    Position_UM = open_file("UM_PositionsInGenome")
    Dico_UM = dic_position(Position_UM, "UM")
    Info_UM = open_file("UM_New_InfoFile")
    ListeUM = update_info(Info_UM, Dico_UM)
    Position_YE = open_file("YE_PositionsInGenome")
    Dico_YE = dic_position(Position_YE, "YE")
    Info_YE = open_file("YE_New_InfoFile")
    ListeYE = update_info(Info_YE, Dico_YE)
    Position_Zamb = open_file("Zamb_PositionsInGenome")
    Dico_Zamb = dic_position(Position_Zamb, "Zamb")
    Info_Zamb = open_file("Zamb_New_InfoFile")
    ListeZamb = update_info(Info_Zamb, Dico_Zamb)
    build_final_file(ListeAK5, ListeDK5, ListeGI5, ListeSW5, ListeUM, ListeYE, ListeZamb, "AllDeNovo_Info")
















