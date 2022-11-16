import cv2
import os
os.chdir("path to folder")
import random


""" docstring

    This program concatenate in one single alignemnts the around 13000 aligned orthogroups
    
    input
    ------------
    - All alignements
    
    output
    ------------
    - Concatenated alignment
    
"""


# GetAllFile Names
def open_file(NameFile):
    F=open(NameFile, "r")
    L=F.readlines()
    return L  


def make_list_names(File):
    # This function builds a list of all proteins names
    L = []
    for i in File:
        Name = i.split("\n")[0]
        NewName = "Ali_"+Name
        L.append(NewName)
    return L


def build_concatd_alignment(ListeNamesAli):
    # This function builds lists of concatenated alignments in each line
    ListeAK5 = []
    ListeDK5 = []
    ListeGI5 = []
    ListeSW5 = []
    ListeUM = []
    ListeYE = []
    ListeZamb = []
    for Name in ListeNamesAli:
        F = open_file(Name)
        NamePop = ""
        Seq = ""
        for i in F:
            ligne = i.split("\n")[0]
            if ligne[0] == ">":
                if NamePop != "":
                    if NamePop == "AK5":
                        ListeAK5.append(Seq)
                    elif NamePop == "DK5":
                        ListeDK5.append(Seq)
                    elif NamePop == "GI5":
                        ListeGI5.append(Seq)
                    elif NamePop == "SW5":
                        ListeSW5.append(Seq)
                    elif NamePop == "UM":
                        ListeUM.append(Seq)
                    elif NamePop == "YE":
                        ListeYE.append(Seq)
                    elif NamePop == "Zamb":
                        ListeZamb.append(Seq)
                    Seq = ""
                NamePop = ligne.split("__")[1]
            else:
                Seq+=ligne
                #print ligne
        if NamePop == "AK5":
            ListeAK5.append(Seq)
        elif NamePop == "DK5":
            ListeDK5.append(Seq)
        elif NamePop == "GI5":
            ListeGI5.append(Seq)
        elif NamePop == "SW5":
            ListeSW5.append(Seq)
        elif NamePop == "UM":
            ListeUM.append(Seq)
        elif NamePop == "YE":
            ListeYE.append(Seq)
        elif NamePop == "Zamb":
            ListeZamb.append(Seq)
        Seq = ""
        NamePop = ""
    return ListeAK5, ListeDK5, ListeGI5, ListeSW5, ListeUM, ListeYE, ListeZamb


def concat(ListeAK5, ListeDK5, ListeGI5, ListeSW5, ListeUM, ListeYE, ListeZamb):
    # This function concatenates the aligned proteins
    SeqAK5 = ""
    SeqDK5 = ""
    SeqGI5 = ""
    SeqSW5 = ""
    SeqUM = ""
    SeqYE = ""
    SeqZamb = ""
    for i in ListeAK5:
        SeqAK5+=i
    for i in ListeDK5:
        SeqDK5+=i
    for i in ListeGI5:
        SeqGI5+=i
    for i in ListeSW5:
        SeqSW5+=i
    for i in ListeUM:
        SeqUM+=i
    for i in ListeYE:
        SeqYE+=i
    for i in ListeZamb:
        SeqZamb+=i
    F = open("PopProtConcatAli.fa", "w")
    F.write(">AK5"+"\n")
    F.write(SeqAK5+"\n")
    F.write(">DK5"+"\n")
    F.write(SeqDK5+"\n")
    F.write(">GI5"+"\n")
    F.write(SeqGI5+"\n")
    F.write(">SW5"+"\n")
    F.write(SeqSW5+"\n")
    F.write(">UM"+"\n")
    F.write(SeqUM+"\n")
    F.write(">YE"+"\n")
    F.write(SeqYE+"\n")
    F.write(">Zamb"+"\n")
    F.write(SeqZamb+"\n")
    F.close()


def main_function():
    FileName = open_file("AllFileNames.txt")
    ListeNamesAli = make_list_names(FileName)
    ListeAK5, ListeDK5, ListeGI5, ListeSW5, ListeUM, ListeYE, ListeZamb = build_concatd_alignment(ListeNamesAli)
    concat(ListeAK5, ListeDK5, ListeGI5, ListeSW5, ListeUM, ListeYE, ListeZamb)












