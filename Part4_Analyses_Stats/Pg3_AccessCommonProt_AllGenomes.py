import cv2
import os
os.chdir("path to folder")
import random


""" docstring

    This program retrieve all proteins in commun between genomes and rename them
    
    input
    ------------
    - protein lists of 7 lines
    
    output
    ------------
    - Renamed common proteins
    
"""


# Build the list with corresponding proteins
def open_file(NameFile):
    F = open(NameFile, "r")
    L = F.readlines()
    return L  


def access_protein_name(File):
    ListeLongName = []
    ListeShortName = []
    for i in File: 
        LongName = i.split(" ")[0]
        if LongName[0] == ">":
            LongName = LongName[1:]
            ShortName = LongName.split(".")[0]
            ListeLongName.append(LongName)
            ListeShortName.append(ShortName)
    return ListeLongName,ListeShortName


def get_commun_names(L1,L2,L3,L4,L5,L6,L7):
    ListeNamesCommun = []
    for i in L1:
        if i in L2 and i in L3 and i in L4 and i in L5 and i in L6 and i in L7:
            ListeNamesCommun.append(i)
    print len(ListeNamesCommun)
    return ListeNamesCommun


# Make All files for alignments
def retreive_proteins_and_rename_them(File):
    Dico = {}
    GeneName = ""
    Prot = ""
    for i in File:
        ligne = i.split("\n")[0]
        if ligne[0] == ">":
            if GeneName!="":
                Dico[GeneName] = Prot[0:len(Prot)-1]
                Prot = ""
            LongName = ligne.split(" ")[0]
            LongName = LongName[1:]
            GeneName = LongName
        else:
            Prot+=ligne
    Dico[GeneName] = Prot[0:len(Prot)-1]
    return Dico


def build_final_file(DicoProtAK5, DicoProtDK5, DicoProtGI5, DicoProtSW5, DicoProtUM, DicoProtYE, DicoProtZamb, ListeCommunLongNames):
    ListeNameFiles = []
    for i in ListeCommunLongNames:
        NameFile = "TheProt_"+i+".fa"
        F = open(NameFile, "w")
        F.write(">"+i+"__"+"AK5"+"\n")
        F.write(DicoProtAK5[i]+"\n")
        F.write(">"+i+"__"+"DK5"+"\n")
        F.write(DicoProtDK5[i]+"\n")
        F.write(">"+i+"__"+"GI5"+"\n")
        F.write(DicoProtGI5[i]+"\n")
        F.write(">"+i+"__"+"SW5"+"\n")
        F.write(DicoProtSW5[i]+"\n")
        F.write(">"+i+"__"+"UM"+"\n")
        F.write(DicoProtUM[i]+"\n")
        F.write(">"+i+"__"+"YE"+"\n")
        F.write(DicoProtYE[i]+"\n")
        F.write(">"+i+"__"+"Zamb"+"\n")
        F.write(DicoProtZamb[i]+"\n")
        F.close()
        ListeNameFiles.append(NameFile)
    F = open("AllFileNames.txt", "w")
    for i in ListeNameFiles:
        F.write(i+"\n")
    F.close()
            

def main_function():
    DataAK5 = open_file("ProteinsAK5.fa")
    ListeLongNameAK5,ListeShortNameAK5 = access_protein_name(DataAK5)
    DataDK5 = open_file("ProteinsDK5.fa")
    ListeLongNameDK5,ListeShortNameDK5 = access_protein_name(DataDK5)
    DataGI5 = open_file("ProteinsGI5.fa")
    ListeLongNameGI5,ListeShortNameGI5 = access_protein_name(DataGI5)
    DataSW5 = open_file("ProteinsSW5.fa")
    ListeLongNameSW5,ListeShortNameSW5 = access_protein_name(DataSW5)
    DataUM = open_file("ProteinsUM.fa")
    ListeLongNameUM,ListeShortNameUM = access_protein_name(DataUM)
    DataYE = open_file("ProteinsYE.fa")
    ListeLongNameYE,ListeShortNameYE = access_protein_name(DataYE)
    DataZamb = open_file("ProteinsZamb.fa")
    ListeLongNameZamb,ListeShortNameZamb = access_protein_name(DataZamb)
    ListeCommunLongNames = get_commun_names(ListeLongNameAK5,ListeLongNameDK5,ListeLongNameGI5,ListeLongNameSW5,ListeLongNameUM,ListeLongNameYE,ListeLongNameZamb)
    get_commun_names(ListeShortNameAK5,ListeShortNameDK5,ListeShortNameGI5,ListeShortNameSW5,ListeShortNameUM,ListeShortNameYE,ListeShortNameZamb)
    DicoProtAK5 = retreive_proteins_and_rename_them(DataAK5)
    DicoProtDK5 = retreive_proteins_and_rename_them(DataDK5)
    DicoProtGI5 = retreive_proteins_and_rename_them(DataGI5)
    DicoProtSW5 = retreive_proteins_and_rename_them(DataSW5)
    DicoProtUM = retreive_proteins_and_rename_them(DataUM)
    DicoProtYE = retreive_proteins_and_rename_them(DataYE)
    DicoProtZamb = retreive_proteins_and_rename_them(DataZamb)
    build_final_file(DicoProtAK5, DicoProtDK5, DicoProtGI5, DicoProtSW5, DicoProtUM, DicoProtYE, DicoProtZamb, ListeCommunLongNames)





