import cv2
import os
os.chdir("path to folder")
import random


""" This program builds orthogroup results according to Definition 3

	input
	------------
	- Info file with unspliced positions of all de novo transcripts
	
	output
	------------
	- Orthogroups according to Definition 3
"""


def open_file(NameFile):
    F=open(NameFile, "r")
    L=F.readlines()
    return L  


def build_dic_de_novo(File):
	L = []
	D = {}
	for i in File[1:]:
		ligne = i.split(",")
		Name = ligne[6]
		Chrom = ligne[2]
		Start = int(ligne[3])
		Stop = int(ligne[4])
		Signe = ligne[5]
		L.append(Name)
		NewL = [Chrom, Start, Stop, Signe]
		D[Name] = NewL
	return L, D
	

def build_orthogroups(LdeNovo, DdeNovo, Interval):
	DicoOrthogroups = {}
	Compteur = 0
	while len(LdeNovo)>0:
		Compteur += 1
		LeaderName = LdeNovo[0]
		LeaderCharacteristics = DdeNovo[LeaderName]
		NewListe = []
		NewListe.append(LeaderName)
		for i in LdeNovo[1:]:
			CandidateName = i
			CandidateCharacteristics = DdeNovo[CandidateName]
			if LeaderCharacteristics[0] == CandidateCharacteristics[0]:
				if LeaderCharacteristics[3] == CandidateCharacteristics[3] or LeaderCharacteristics[3] == "." or CandidateCharacteristics[3] == ".":
					if abs(CandidateCharacteristics[1]-LeaderCharacteristics[1]) <= Interval:
						if abs(CandidateCharacteristics[2]-LeaderCharacteristics[2]) <= Interval:
							NewListe.append(CandidateName)
		for j in NewListe:
			LdeNovo.remove(j)
		Name = "Orthogroups"+str(Compteur)
		DicoOrthogroups[Name] = NewListe
	print len(DicoOrthogroups)
	return DicoOrthogroups


def build_final_file(Name, Dico):
	F = open(Name, "w")
	for i in Dico.keys():
		Name = i
		Liste = Dico[i]
		Size = len(Liste)
		F.write(Name+","+str(Size))
		for j in Liste:
			F.write(","+j)
		F.write("\n")
	F.close()


def main_function():
	AllDeNovo = open_file("AllDeNovo_Info")
	LdeNovo, DdeNovo = build_dic_de_novo(AllDeNovo)
	DicoOrthogroups = build_orthogroups(LdeNovo, DdeNovo, 500)
	build_final_file("DeNovoTranscriptsGroups_Rst3", DicoOrthogroups)


