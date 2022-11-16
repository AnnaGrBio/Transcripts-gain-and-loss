import cv2
import os
os.chdir("path to folder")
import random


""" This progra format orthogroups from Definition 1, 2 and 3, according to the genomic position of the DeNovoTranscriptsGroups_Rst1

	input
	------------
	- Info file all de novo transcripts
	- Orthogroups from Definition 1, 2 and 3
	
	output
	------------
	- Orthogroups from Definition 1, 2 and 3 formated.
	
"""


def open_file(NameFile):
    F=open(NameFile, "r")
    L=F.readlines()
    return L  


def MakeDico(File, Info):
	D = {}
	for i in File:
		ligne = i.split(",")
		NameOrthogroup = ligne[0]
		Liste = []
		for j in ligne[2:]:
			NamePop = j.split("_")[0]
			if NamePop not in Liste:
				Liste.append(NamePop)
		NameRef = ligne[2].split("\n")
		NameRef = NameRef[0]
		for k in Info:
			lala = k.split(",")
			if lala[6] == NameRef:
				Cat = lala[14].split("\n")[0]
				Liste.append(Cat)
				break
		D[NameOrthogroup]=Liste
	return D


def MakeListeToSort():
	L = []
	C = 1
	for i in range(50000):
		Name = "Orthogroups"+str(C)
		L.append(Name)
		C +=1
	return L


def MakeFinalFile(Dico, Liste, Name):
	F = open(Name, "w")
	for i in Liste:
		if i in Dico.keys():
			SubListe = Dico[i]
			F.write(i)
			for j in SubListe:
				F.write(","+j)
			F.write("\n")
		else:
			break
	F.close()


def main_function():
	SortingList = MakeListeToSort()
	Info = open_file("AllDeNovo_Info")
	FileRst2 = open_file("DeNovoTranscriptsGroups_Rst1")
	DicoRst2 = MakeDico(FileRst2, Info)
	MakeFinalFile(DicoRst2, SortingList, "DeNovoTranscript_Rst1_FormatedWithCat.txt")
	FileRst3 = open_file("DeNovoTranscriptsGroups_Rst2")
	DicoRst3 = MakeDico(FileRst3, Info)
	MakeFinalFile(DicoRst3, SortingList, "DeNovoTranscript_Rst2_FormatedWithCat.txt")
	FileRst4 = open_file("DeNovoTranscriptsGroups_Rst3")
	DicoRst4 = MakeDico(FileRst4, Info)
	MakeFinalFile(DicoRst4, SortingList, "DeNovoTranscript_Rst3_FormatedWithCat.txt")


