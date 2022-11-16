import cv2
import os
from Bio.Seq import Seq
os.chdir("path to folder")
import random


""" docstring

	This program reads into BLAST output and retrieves all transcripts with no hits
	
	input
	------------
	BLAST output
	
	output
	------------
	name of BLASTs output sequences without hits
"""


def open_file(NameFile):
    F=open(NameFile, "r")
    L=F.readlines()
    return L  


def find_no_hit(F, Evalue):
	# This function find transcripts without hits and store their name in a list
	NbTotalQuery = 0
	ListeFinale = []
	for i in range(0,len(F)):
		if F[i][0:6] == "Query=":
			NbTotalQuery+=1
			#print NbTotalQuery
			if F[i+5] == "***** No hits found *****"+"\n":
				Nom =  F[i][7:]
				ListeFinale.append(Nom)
			else:
				#print F[i+6]
				lala = F[i+6]
				lili = lala.split()
				StrValue = lili[len(lili)-1]
				FValue = float(StrValue)
				if float(FValue)>float(Evalue):
					Nom =  F[i][7:]
					ListeFinale.append(Nom)
	print ("Nb initial ORF : "+str(NbTotalQuery))
	print ("Nb Total de novo : "+str(len(ListeFinale)))
	return ListeFinale


def build_final_file(Name, Liste):
	# This function builds the final file of name of transcripts witout hits
	F = open(Name, "w")
	for i in Liste:
		F.write(i)
	F.close()


def main_function():
	# Main function, input have to be modified for each line
	F = open_file("AK5_Trancript_Forward_NonDipteran")
	print ("******************")
	print ("FileOpened")
	print ("******************")
	Evalue = 0.01
	Liste = find_no_hit(F, Evalue)
	build_final_file("AK5_DeNovo_Trancript_Forward_NonDipteran", Liste)
