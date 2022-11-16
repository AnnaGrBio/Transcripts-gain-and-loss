import cv2
import os
os.chdir("path to folder")
import random


""" docstring

	This program does formatting of orthogroups from definition 1, 2 and 3
	
	input
	------------
	- Orthogroups from definition 1, 2 and 3
	
	output
	------------
	- Orthogropus formated for Definition 1, 2 and 3
"""


def open_file(NameFile):
    F = open(NameFile, "r")
    L = F.readlines()
    return L  


def build_dictionary(File):
	D = {}
	for i in File:
		ligne = i.split(",")
		NameOrthogroup = ligne[0]
		Liste = []
		for j in ligne[2:]:
			NamePop = j.split("_")[0]
			if NamePop not in Liste:
				Liste.append(NamePop)
		D[NameOrthogroup]=Liste
	return D


def build_list_to_sort():
	L = []
	C = 1
	for i in range(50000):
		Name = "Orthogroups"+str(C)
		L.append(Name)
		C +=1
	return L


def build_final_file(Dico, Liste, Name):
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
	SortingList = build_list_to_sort()
	FileRst2 = open_file("DeNovoTranscriptsGroups_Rst1")
	DicoRst2 = build_dictionary(FileRst2)
	build_final_file(DicoRst2, SortingList, "DeNovoTranscriptRst1Formated")
	FileRst3 = open_file("DeNovoTranscriptsGroups_Rst2")
	DicoRst3 = build_dictionary(FileRst3)
	build_final_file(DicoRst3, SortingList, "DeNovoTranscriptRst2Formated")
	FileRst4 = open_file("DeNovoTranscriptsGroups_Rst3")
	DicoRst4 = build_dictionary(FileRst4)
	build_final_file(DicoRst4, SortingList, "DeNovoTranscriptRst3Formated")


