import cv2
import os
os.chdir("path to folder")
import random


""" docstring

    This program build the first draft of transcripts orthogroups basd on the BLAST results
    
    input
    ------------
    - BLAST results from BLAST between the 7 lines
    
    output
    ------------
    - Orthogroups draft 1
    
"""


def open_file(NameFile):
    F=open(NameFile, "r")
    L=F.readlines()
    return L  


def build_orthogroups(File):
    # This function builds the orthogroups based on BLAST hits
    DicoGroups = {}
    c = 1
    for i in File:
        ligne = i.split("\n")[0]
        ligne = ligne.split(",")
        if len(DicoGroups) == 0:
            Groupe = "Group"+str(c)
            c+=1
            DicoGroups[Groupe] = [ligne[0], ligne[1]]
        else:
            Binome1 = ligne[0]
            Binome2 = ligne[1]
            if Binome2 == "NotHit":
                Binome2 = Binome1
            Attributed = False
            for NameGroupes in DicoGroups.keys():
                MemberOfGroup = DicoGroups[NameGroupes]
                if Binome1 in MemberOfGroup and Binome2 in MemberOfGroup:

                    Attributed = True
                    break
                elif Binome1 in MemberOfGroup and Binome2 not in MemberOfGroup:

                    DicoGroups[NameGroupes].append(Binome2)
                    Attributed = True
                    break
                elif Binome2 in MemberOfGroup and Binome1 not in MemberOfGroup:

                    DicoGroups[NameGroupes].append(Binome1)
                    Attributed = True
                    break
            if Attributed == False:
                Groupe = "Group"+str(c)
                c+=1
                if Binome1 == Binome2:
                    DicoGroups[Groupe] = [Binome1]
                else:

                    DicoGroups[Groupe] = [Binome1, Binome2]
    return DicoGroups


def merge_dics(Dico):
    # This function merges the dictionaries of orthogropus transcripts
    print (len(Dico))
    for NameGroupes in Dico.keys():
        if NameGroupes in Dico.keys():
            Liste = Dico[NameGroupes]
            for Name in Liste:
                for NewNames in Dico.keys():
                    if NewNames != NameGroupes:
                        NewListe = Dico[NewNames]
                        if Name in NewListe:
                            for k in NewListe:
                                if k not in Liste:
                                    Dico[NameGroupes].append(k)
                            Dico.pop(NewNames)
                        
    print len(Dico)
    return Dico
        
    
def build_orthogroup_file(D, NameFile):
    # This function builds the file of orthogroups
    print ("The total number of groups is: "+str(len(D)))
    F = open(NameFile, "w")
    Bigest1 = 0
    Bigest2 = 0
    Bigest3 = 0
    for i in D.keys():
        F.write(i+",")
        Liste = D[i]
        F.write(str(len(Liste)))
        for j in Liste:
            F.write(","+j)
        F.write("\n")
        if len(Liste)>Bigest1:
            Bigest1 = len(Liste)
        elif len(Liste)>Bigest2 and len(Liste)<Bigest1:
            Bigest2 = len(Liste)
        elif len(Liste)>Bigest3 and len(Liste)<Bigest2:
            Bigest3 = len(Liste)
    print (Bigest1)
    print (Bigest2)
    print (Bigest3)
    F.close()
        

def main_function():
    File = open_file("AllBLAST.rst")
    D = build_orthogroups(File)
    D = merge_dics(D)
    D = merge_dics(D)
    build_orthogroup_file(D, "Allgroups.rst")
