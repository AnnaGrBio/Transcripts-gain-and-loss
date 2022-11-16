import cv2
import os
os.chdir("path to folder")
import random
from Bio import SeqIO
from Bio import motifs
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
import random
from multiprocessing import Pool


""" docstring

    This program builds orthogroups according to Definition 4
    
    input
    ------------
    - All de novo info 
    
    output
    ------------
    - Orthogroups acording to Definition 4
"""


def open_file(NameFile):
    F=open(NameFile, "r")
    L=F.readlines()
    return L  


# Step 1 distribute by population and orientation
def get_position_unspliced_transcripts(NamePop, File):
    # This function access and store the position in genomes of unspliced transcripts
    ListeForward = []
    ListeReverse = []
    for i in File[1:]:
        ligne = i.split(",")
        last = ligne[14].split("\n")[0]
        if last !="ExonInside":
            TranscriptName = ligne[0]
            Chromosome = ligne[2]
            Start = int(ligne[3])
            End = int(ligne[4])
            Orientation = ligne[5]
            Pop = TranscriptName.split("_")[0]
            if Pop == NamePop:
                Pos = [Chromosome, Start, End]
                if Orientation == "+" or Orientation == ".":
                    ListeForward.append(Pos)
                else:
                    ListeReverse.append(Pos)
    return ListeForward, ListeReverse


# Step 2 Merge Overlap between traanscripts
def divide_list_into_sublists(Liste):
    # This function split transcripts according to the chromosome to which they belong
    Data2L = []
    Data2R = []
    Data3L = []
    Data3R = []
    Data4 = []
    DataMito = []
    DataX = []
    DataY = []
    for i in Liste:
        NewListe = [i[1],i[2]]
        if i[0] == "2L_Chromosome":
            Data2L.append(NewListe)
        elif i[0] == "2R_Chromosome":
            Data2R.append(NewListe)
        elif i[0] == "3L_Chromosome":
            Data3L.append(NewListe)        
        elif i[0] == "3R_Chromosome":
            Data3R.append(NewListe)
        elif i[0] == "4_Chromosome":
            Data4.append(NewListe)
        elif i[0] == "X_Chromosome":
            DataX.append(NewListe)
        elif i[0] == "Y_Chromosome":
            DataY.append(NewListe)
        elif i[0] == "mitochondrion_Chromosome":
            DataMito.append(NewListe)
    return [Data2L,Data2R,Data3L,Data3R,Data4,DataMito,DataX,DataY]


def sort_list_a(Liste):
    # This function sort a list with overlap in increasing order
    NewList = []
    for i in Liste:
        if len(NewList) == 0:
            NewList.append(i)
        else:
            Placement = False
            Start_i = i[0]
            End_i = i[1]
            for j in range(len(NewList)):
                ValueStart = NewList[j][0]
                ValueEnd = NewList[j][1]
                if Start_i>ValueEnd:
                    continue
                elif Start_i<ValueStart and End_i<ValueStart:
                    NewList.insert(j,i)
                    Placement = True
                    break
                elif Start_i<=ValueStart and End_i>=ValueStart and End_i<=ValueEnd:
                    NewList[j][0] = Start_i
                    Placement = True
                    break
                elif Start_i>=ValueStart and Start_i<=ValueEnd and End_i>=ValueEnd:
                    NewList[j][1] = End_i
                    Placement = True
                    break
                elif Start_i<ValueStart and End_i>ValueEnd:
                    NewList[j][0] = Start_i
                    NewList[j][1] = End_i
                    Placement = True
                    break
                elif Start_i>ValueStart and End_i<ValueEnd:
                    Placement = True
                    break
            if Placement == False:
                NewList.append(i)
    print (len(NewList))
    return NewList
                        
            
def sort_liste_b(Liste):
    # This function creates a dictionary with transcripts per chromosome in increasing order of position in the chromosome
    ListsDivided = divide_list_into_sublists(Liste)
    NewListDivided = []
    for i in ListsDivided:
        if len(i)>0:
            print len(i)
            T1 = sort_list_a(i) #1. Sort the list in increasing order and revove overlaps
            T2 = sort_list_a(T1)
            T3 = sort_list_a(T2)
            NewListDivided.append(T3)
            print ("*********************")
        else:
            NewListDivided.append(i)
            
    DicoSorted = {"2L":NewListDivided[0], "2R":NewListDivided[1],"3L":NewListDivided[2],"3R":NewListDivided[3],"4":NewListDivided[4],"Mito":NewListDivided[5],"X":NewListDivided[6],"Y":NewListDivided[7]}
    return DicoSorted
    

# Step 3 compare coverage

SizeWholeChrom2L = 23513712
SizeWholeChrom2R = 25286936
SizeWholeChrom3L = 28110227
SizeWholeChrom3R = 32079331
SizeWholeChrom4 = 1348131
SizeWholeChromMito = 19524
SizeWholeChromX = 23542271
SizeWholeChromY = 3667352
SizeTotal = 137567484
DicoCorresChrom = {"2L":SizeWholeChrom2L, "2R":SizeWholeChrom2R,"3L":SizeWholeChrom3L,"3R":SizeWholeChrom3R,"4":SizeWholeChrom4,"Mito":SizeWholeChromMito,"X":SizeWholeChromX,"Y":SizeWholeChromY}


def make_list_all_poss(Liste):
    # This function makes the list of all positions
    NewListe = []
    for i in Liste:
        for number in range(i[0],i[1]+1):
            NewListe.append(number)
    return NewListe


def implement_dic(DicoChromSpecific,chrom):
    NewDico = {}
    while len(DicoChromSpecific)>0:
        Combinaison = DicoChromSpecific[DicoChromSpecific.keys()[0]]
        Compteur = 0
        for i in DicoChromSpecific.keys():
            if DicoChromSpecific[i] == Combinaison:
                Compteur += 1
                del DicoChromSpecific[i]
        NewDico[Combinaison] = Compteur
    return (NewDico)
                

def build_dic_orthogroups(DicoAK5,DicoDK5,DicoGI5,DicoSW5,DicoUM,DicoYE,DicoZamb,ListeChrom):
    # This function builds dictionary of orthogroups
    ListeFinalDico = []
    FinalDicoPerChrom = {}
    for chrom in ListeChrom:
        DicoChromSpecific = {}
        ListeTranscAK5=DicoAK5[chrom]
        PosAK5 = make_list_all_poss(ListeTranscAK5)
        ListeTranscDK5=DicoDK5[chrom]
        PosDK5 = make_list_all_poss(ListeTranscDK5)
        ListeTranscGI5=DicoGI5[chrom]
        PosGI5 = make_list_all_poss(ListeTranscGI5)
        ListeTranscSW5=DicoSW5[chrom]
        PosSW5 = make_list_all_poss(ListeTranscSW5)
        ListeTranscUM=DicoUM[chrom]
        PosUM = make_list_all_poss(ListeTranscUM)
        ListeTranscYE=DicoYE[chrom]
        PosYE = make_list_all_poss(ListeTranscYE)
        ListeTranscZamb=DicoZamb[chrom]
        PosZamb = make_list_all_poss(ListeTranscZamb)
        print ("List all pos done")
        print ("Start with AK5")
        print ("Size AK5 : "+str(len(PosAK5)))
        c = 0
        for Position in PosAK5:
            #c+=1
            #print c
            SetPops="AK5"
            if Position in PosDK5:
                SetPops+=",DK5"
                PosDK5.remove(Position)
            if Position in PosGI5:
                SetPops+=",GI5"
                PosGI5.remove(Position)
            if Position in PosSW5:
                SetPops+=",SW5"
                PosSW5.remove(Position)
            if Position in PosUM:
                SetPops+=",UM"
                PosUM.remove(Position)
            if Position in PosYE:
                SetPops+=",YE"
                PosYE.remove(Position)
            if Position in PosZamb:
                SetPops+=",Zamb"
                PosZamb.remove(Position)
            DicoChromSpecific[Position] = SetPops
        print ("Start with DK5")
        print ("Size DK5 : "+str(len(PosAK5)))
        c = 0
        for Position in PosDK5:
            #c+=1
            #print c
            SetPops="DK5"
            if Position in PosGI5:
                SetPops+=",GI5"
                PosGI5.remove(Position)
            if Position in PosSW5:
                SetPops+=",SW5"
                PosSW5.remove(Position)
            if Position in PosUM:
                SetPops+=",UM"
                PosUM.remove(Position)
            if Position in PosYE:
                SetPops+=",YE"
                PosYE.remove(Position)
            if Position in PosZamb:
                SetPops+=",Zamb"
                PosZamb.remove(Position)
            DicoChromSpecific[Position] = SetPops 
        for Position in PosGI5:
            SetPops="GI5"
            if Position in PosSW5:
                SetPops+=",SW5"
                PosSW5.remove(Position)
            if Position in PosUM:
                SetPops+=",UM"
                PosUM.remove(Position)
            if Position in PosYE:
                SetPops+=",YE"
                PosYE.remove(Position)
            if Position in PosZamb:
                SetPops+=",Zamb"
                PosZamb.remove(Position)
            DicoChromSpecific[Position] = SetPops
        for Position in PosSW5:
            SetPops="SW5"
            if Position in PosUM:
                SetPops+=",UM"
                PosUM.remove(Position)
            if Position in PosYE:
                SetPops+=",YE"
                PosYE.remove(Position)
            if Position in PosZamb:
                SetPops+=",Zamb"
                PosZamb.remove(Position)
            DicoChromSpecific[Position] = SetPops
        for Position in PosUM:
            SetPops="UM"
            if Position in PosYE:
                SetPops+=",YE"
                PosYE.remove(Position)
            if Position in PosZamb:
                SetPops+=",Zamb"
                PosZamb.remove(Position)
            DicoChromSpecific[Position] = SetPops
        for Position in PosYE:
            SetPops="YE"
            if Position in PosZamb:
                SetPops+=",Zamb"
                PosZamb.remove(Position)
            DicoChromSpecific[Position] = SetPops
        for Position in PosZamb:
            SetPops="Zamb"
            DicoChromSpecific[Position] = SetPops
        FinalDico = implement_dic(DicoChromSpecific,chrom)
    return FinalDico
    
        
def merge_all_informations(ListeFinalDicoForward,ListeFinalDicoReverse):
    # THis function merges information of all chromosomes
    SizeTotal = 137567484
    ListeAllPossibleCombis = []
    for dico in ListeFinalDicoForward:
        for combi in dico.keys():
            if combi not in ListeAllPossibleCombis:
                ListeAllPossibleCombis.append(combi)
    FinalTotalDicoNumbers = {}
    for dico in ListeFinalDicoReverse:
        for combi in dico.keys():
            if combi not in ListeAllPossibleCombis:
                ListeAllPossibleCombis.append(combi)
    for combi in ListeAllPossibleCombis:
        compteurTotal = 0
        for dico in ListeFinalDicoForward:
            if combi in dico.keys():
                Nb = dico[combi]
                compteurTotal += Nb
        for dico in ListeFinalDicoReverse:
            if combi in dico.keys():
                Nb = dico[combi]
                compteurTotal += Nb
        FinalTotalDicoNumbers[combi] = compteurTotal
    # Get percentage of occurence of each combi
    FinalTotalDicoPercentageRecouvrament = {}
    for combi in FinalTotalDicoNumbers.keys():
        Number = FinalTotalDicoNumbers[combi]
        Percentage = float(100)*float(Number)/float(SizeTotal*2)
        FinalTotalDicoPercentageRecouvrament[combi] = Percentage
    return FinalTotalDicoNumbers, FinalTotalDicoPercentageRecouvrament
    
    
def wrapper(x):
	DicoAK5,DicoDK5,DicoGI5,DicoSW5,DicoUM,DicoYE,DicoZamb,ListeChrom = x
	D = build_dic_orthogroups(DicoAK5,DicoDK5,DicoGI5,DicoSW5,DicoUM,DicoYE,DicoZamb,ListeChrom)
	return D


def Multiproc(L):
	Taille = len(L)
	p = Pool(Taille)
	ListeFinalDicoForward = p.map(wrapper, L)
	p.close()
	p.join()
	print len(ListeFinalDicoForward)
	return ListeFinalDicoForward
	#join the 60 dic from the Results list


def main_function():
    print ("Start Step 1 : store de novo transcripts by population and orientation")
    BigFile = open_file("AllDeNovo_Info")
    ListeForwardAK5, ListeReverseAK5 = get_position_unspliced_transcripts("AK5", BigFile)
    print(len(ListeForwardAK5))
    print(len(ListeReverseAK5))
    ListeForwardDK5, ListeReverseDK5 = get_position_unspliced_transcripts("DK5", BigFile)
    print(len(ListeForwardDK5))
    print(len(ListeReverseDK5))
    ListeForwardGI5, ListeReverseGI5 = get_position_unspliced_transcripts("GI5", BigFile)
    print(len(ListeForwardGI5))
    print(len(ListeReverseGI5))
    ListeForwardSW5, ListeReverseSW5 = get_position_unspliced_transcripts("SW5", BigFile)
    print(len(ListeForwardSW5))
    print(len(ListeReverseSW5))
    ListeForwardUM, ListeReverseUM = get_position_unspliced_transcripts("UM", BigFile)
    print(len(ListeForwardUM))
    print(len(ListeReverseUM))
    ListeForwardYE, ListeReverseYE = get_position_unspliced_transcripts("YE", BigFile)
    print(len(ListeForwardYE))
    print(len(ListeReverseYE))
    ListeForwardZamb, ListeReverseZamb = get_position_unspliced_transcripts("Zamb", BigFile)
    print(len(ListeForwardZamb))
    print(len(ListeReverseZamb))
    print ("Start Step 2 : Merge Overlap between transcripts in each population by orientation")
    DicoAK5forward = sort_liste_b(ListeForwardAK5)
    DicoAK5reverse = sort_liste_b(ListeReverseAK5)
    DicoDK5forward = sort_liste_b(ListeForwardDK5)
    DicoDK5reverse = sort_liste_b(ListeReverseDK5)
    DicoGI5forward = sort_liste_b(ListeForwardGI5)
    DicoGI5reverse = sort_liste_b(ListeReverseGI5)
    DicoSW5forward = sort_liste_b(ListeForwardSW5)
    DicoSW5reverse = sort_liste_b(ListeReverseSW5)
    DicoUMforward = sort_liste_b(ListeForwardUM)
    DicoUMreverse = sort_liste_b(ListeReverseUM)
    DicoYEforward = sort_liste_b(ListeForwardYE)
    DicoYEreverse = sort_liste_b(ListeReverseYE)
    DicoZambforward = sort_liste_b(ListeForwardZamb)
    DicoZambreverse = sort_liste_b(ListeReverseZamb)
    print ("Start Step 3 : Make orthogroup and count number of nucleotide in common per orthogroup, in each chrom and orientation")
    ListeChrom1 = ["2L"]
    ListeChrom2 = ["2R"]
    ListeChrom3 = ["3L"]
    ListeChrom4 = ["3R"]
    ListeChrom5 = ["4"]
    ListeChrom6 = ["Mito"]
    ListeChrom7 = ["X"]
    ListeChrom8 = ["Y"]
    ListeToParse1 = [[DicoAK5forward,DicoDK5forward,DicoGI5forward,DicoSW5forward,DicoUMforward,DicoYEforward,DicoZambforward,ListeChrom1],[DicoAK5forward,DicoDK5forward,DicoGI5forward,DicoSW5forward,DicoUMforward,DicoYEforward,DicoZambforward,ListeChrom2],[DicoAK5forward,DicoDK5forward,DicoGI5forward,DicoSW5forward,DicoUMforward,DicoYEforward,DicoZambforward,ListeChrom3],[DicoAK5forward,DicoDK5forward,DicoGI5forward,DicoSW5forward,DicoUMforward,DicoYEforward,DicoZambforward,ListeChrom4],[DicoAK5forward,DicoDK5forward,DicoGI5forward,DicoSW5forward,DicoUMforward,DicoYEforward,DicoZambforward,ListeChrom5],[DicoAK5forward,DicoDK5forward,DicoGI5forward,DicoSW5forward,DicoUMforward,DicoYEforward,DicoZambforward,ListeChrom6],[DicoAK5forward,DicoDK5forward,DicoGI5forward,DicoSW5forward,DicoUMforward,DicoYEforward,DicoZambforward,ListeChrom7],[DicoAK5forward,DicoDK5forward,DicoGI5forward,DicoSW5forward,DicoUMforward,DicoYEforward,DicoZambforward,ListeChrom8]]
    ListeToParse2 = [[DicoAK5reverse,DicoDK5reverse,DicoGI5reverse,DicoSW5reverse,DicoUMreverse,DicoYEreverse,DicoZambreverse,ListeChrom1],[DicoAK5reverse,DicoDK5reverse,DicoGI5reverse,DicoSW5reverse,DicoUMreverse,DicoYEreverse,DicoZambreverse,ListeChrom2],[DicoAK5reverse,DicoDK5reverse,DicoGI5reverse,DicoSW5reverse,DicoUMreverse,DicoYEreverse,DicoZambreverse,ListeChrom3],[DicoAK5reverse,DicoDK5reverse,DicoGI5reverse,DicoSW5reverse,DicoUMreverse,DicoYEreverse,DicoZambreverse,ListeChrom4],[DicoAK5reverse,DicoDK5reverse,DicoGI5reverse,DicoSW5reverse,DicoUMreverse,DicoYEreverse,DicoZambreverse,ListeChrom5],[DicoAK5reverse,DicoDK5reverse,DicoGI5reverse,DicoSW5reverse,DicoUMreverse,DicoYEreverse,DicoZambreverse,ListeChrom6],[DicoAK5reverse,DicoDK5reverse,DicoGI5reverse,DicoSW5reverse,DicoUMreverse,DicoYEreverse,DicoZambreverse,ListeChrom7],[DicoAK5reverse,DicoDK5reverse,DicoGI5reverse,DicoSW5reverse,DicoUMreverse,DicoYEreverse,DicoZambreverse,ListeChrom8]]
    ListeFinalDicoForward = Multiproc(ListeToParse1)
    ListeFinalDicoReverse = Multiproc(ListeToParse2)
    print ("Start Step 4: Merge dictionaries for final result")
    FinalTotalDicoNumbers, FinalTotalDicoPercentageRecouvrament = merge_all_informations(ListeFinalDicoForward,ListeFinalDicoReverse)
    F = open("CoverageOrthogroupsNUMBERs.rst4", "w")
    for i in FinalTotalDicoNumbers.keys():
        F.write(i+","+str(FinalTotalDicoNumbers[i])+"\n")
    F.close()
    F = open("CoverageOrthogroupsPERCENTAGEs.rst4", "w")
    for i in FinalTotalDicoPercentageRecouvrament.keys():
        F.write(i+","+str(FinalTotalDicoPercentageRecouvrament[i])+"\n")
    F.close()
