import cv2
import os
os.chdir("path to folder")
import random


""" This program builds transcripts orthogroups according to Definition 1 and Definition 2

    input
    ------------
    - Orthogroup first draft
    - Info file of all de novo transcripts
    
    output
    ------------
    - 1 file with orthogroups according to Definition 1
    - 1 file with orthogroups according to Definition 2
    
"""


def open_file(NameFile):
    F=open(NameFile, "r")
    L=F.readlines()
    return L  


def store_groups(Groups):
    # This function stores groups in a dictionary and returns it
    Dico = {}
    for i in Groups:
        ligne = i.split("\n")[0]
        Name = ligne.split(",")[0]
        Liste = ligne.split(",")
        NewListe = []
        for j in Liste[2:]:
            NewListe.append(j)
        Dico[Name] = NewListe
    return Dico


def assess_transcripts_size(Info):
    # This function assess the size of spliced transcripts
    Min = 10000
    NamMax = ""
    Max = 0
    for i in Info[1:]:
        ligne = i.split(",")
        Beg = int(ligne[3])
        End = int(ligne[4])
        Size = End-Beg
        if Size>Max:
            Max = Size
            NamMax = ligne[6]
        if Size<Min:
            Min = Size
    print ("Max Transcript Size = "+str(Max))
    print NamMax
    print ("Min Transcript Size = "+str(Min))  
    
    
def build_dic_pos(Info):
    # This function estimates the starting point of the transcripts
    Dico = {}
    for i in Info[1:]:
        ligne = i.split(",")
        Name = ligne[6]
        Chrom = ligne[2]
        Beg = int(ligne[3])
        End = int(ligne[4])
        Signe = ligne[5]
        Dico[Name] = [Chrom, Beg, End, Signe]
    return Dico
        
        
def test_if_all_correct(GroupsMembers, DicoPos):
    # This function calculate overlap to make sure the coverage of overlaps are or are not corrects
    Correct = True
    for Transcript in GroupsMembers:
        Data = DicoPos[Transcript]
        Chrom = Data[0]
        Beg = Data[1]
        End = Data[2]
        Overlap = False
        for OtherTransc in GroupsMembers:
            if OtherTransc!=Transcript:
                OtherData = DicoPos[OtherTransc]
                OtherChrom = OtherData[0]
                OtherBeg = OtherData[1]
                OtherEnd = OtherData[2]
                if OtherChrom!=Chrom:
                    Correct = False
                    break
                else:
                    if Beg<=OtherBeg and End>OtherBeg:
                        Overlap = True
                    elif End>=OtherEnd and Beg<OtherEnd:
                        Overlap = True
                    elif Beg>OtherBeg and End<OtherEnd:
                        Overlap = True
                    elif Beg==OtherBeg and End==OtherEnd:
                        Overlap = True
            else:
                Overlap = True
        if Overlap == False:
            Correct = False
        if Correct == False:
            break
    return Correct


def build_new_groups(GroupsMembers,DicoPos):
    # This function starts building new orthogroups based on the previous assessments and informations
    FinalListe = []
    SortedByChrom = []
    Chrom3L = []
    Chrom3R = []
    Chrom2L = []
    Chrom2R = []
    Chrom4 = []
    ChromMito = []
    ChromX = []
    ChromY = []
    for i in GroupsMembers:
        Data = DicoPos[i]
        Chrom = Data[0]
        if Chrom == "2L_Chromosome":
            Chrom2L.append(i)
        elif Chrom == "2R_Chromosome":
            Chrom2R.append(i)
        elif Chrom == "3L_Chromosome":
            Chrom3L.append(i)
        elif Chrom == "3R_Chromosome":
            Chrom3R.append(i)
        elif Chrom == "mitochondrion_Chromosome":
            ChromMito.append(i)
        elif Chrom == "4_Chromosome":
            Chrom4.append(i)
        elif Chrom == "X_Chromosome":
            ChromX.append(i)
        elif Chrom == "Y_Chromosome":
            ChromY.append(i)
    if len(Chrom3L)>0:
        SortedByChrom.append(Chrom3L)
    if len(Chrom3R)>0:
        SortedByChrom.append(Chrom3R)
    if len(Chrom2R)>0:
        SortedByChrom.append(Chrom2R)
    if len(Chrom2L)>0:
        SortedByChrom.append(Chrom2L)
    if len(Chrom4)>0:
        SortedByChrom.append(Chrom4)
    if len(ChromMito)>0:
        SortedByChrom.append(ChromMito)
    if len(ChromX)>0:
        SortedByChrom.append(ChromX)
    if len(ChromY)>0:
        SortedByChrom.append(ChromY)
    NewlySorted = True
    for i in SortedByChrom:
        lala = test_if_all_correct(i, DicoPos)
        if lala == False:
            NewlySorted = False
    if NewlySorted == False:
        print ("Attention!!! Malgres le partage par chromose, il y a un groupe sans overlap")
    return SortedByChrom
        

def rebuild_groups(DicoGroups, DicoPos):
    # This function merge results and re-build correct orthogroups
    NewDico = {}
    c = 1
    NbFalse = 0
    print ("There is initialy a total of "+str(len(DicoGroups)))
    for i in DicoGroups.keys():
        GroupsMembers = DicoGroups[i]
        if len(DicoGroups[i]) == 1:
            NewName = "Orthogroups"+str(c)
            NewDico[NewName] = GroupsMembers
            c+=1
        else:
            #print i
            #print GroupsMembers
            Value = test_if_all_correct(GroupsMembers, DicoPos)
            if Value == True:
                NewName = "Orthogroups"+str(c)
                NewDico[NewName] = GroupsMembers
                c+=1
            else:
                SortedByChrom = build_new_groups(GroupsMembers,DicoPos)
                for Listes in SortedByChrom:
                    NewName = "Orthogroups"+str(c)
                    NewDico[NewName] = Listes
                    c+=1
                #print i
                #print GroupsMembers
                NbFalse+=1
    print ("There is a final number of "+str(len(NewDico)))
    print (str(NbFalse)+" groups were reshuffled")
    return NewDico
                

def split_by_starting_point(GroupsMembers, DicoPos, WindowSize):
    # This function separate overlapping transcripts according to the window of their starting point
    FinalListe = []
    for Transcript in GroupsMembers:
        Data = DicoPos[Transcript]
        Signe = Data[3]
        if Signe == "+" or ".":
            Beg = Data[1]
        elif Signe == "-":
            Beg = Data[2]
        if len(FinalListe) == 0:
            FinalListe.append([Transcript])
        else:
            c = 0
            Placement = False
            for Listes in FinalListe:
                OverlapWithAll = True
                for j in Listes:
                    Data2 = DicoPos[j]
                    Signe2 = Data2[3]
                    if Signe2 == "+" or ".":
                        Beg2 = Data2[1]
                    elif Signe2 == "-":
                        Beg2 = Data2[2]
                    if abs(Beg2-Beg)<=WindowSize:
                        continue
                    else:
                        OverlapWithAll = False
                        break
                if OverlapWithAll == True:
                    FinalListe[c].append(Transcript)
                    Placement = True
                    break
                else:
                    c+=1
            if Placement == False:
                FinalListe.append([Transcript])
    #print FinalListe
    return FinalListe
                
                    
def sort_by_starting_point(DicoFinal1, DicoPos, WindowSize):
    # This function distribute overlapping transcripts according to their starting points
    DicoFinal2 = {}
    c = 1
    for i in DicoFinal1.keys():
        GroupsMembers = DicoFinal1[i]
        if len(DicoFinal1[i]) == 1:
            NewName = "Orthogroups"+str(c)
            DicoFinal2[NewName] = GroupsMembers
            c+=1
        else:
            NewListe = split_by_starting_point(GroupsMembers, DicoPos, WindowSize)
            for Listes in NewListe:
                NewName = "Orthogroups"+str(c)
                DicoFinal2[NewName] = Listes
                c+=1
    print ("There is a final number of "+str(len(DicoFinal2)))
    return DicoFinal2


def calculate_coverage(Beg, End, Beg2, End2, Coverage):
    # This function calculate the percentage of overlap between transcripts
    Cover = False
    SizeT1 = End-Beg    
    SizeT2 = End2-Beg2
    if Beg<=Beg2 and End>=End2:
        Cover = True
    elif Beg<=Beg2 and End>Beg2:
        if End>End2:
            Cover = True
        else:
            CoveredSize = End-Beg2
            PercCover1 = 100*CoveredSize/SizeT1
            PercCover2 = 100*CoveredSize/SizeT2
            if PercCover1>Coverage or PercCover2>Coverage:
                Cover = True
    elif End>=End2 and Beg<End2:
        if Beg<Beg2:
            Cover = True
        else:
            CoveredSize = End2-Beg
            PercCover1 = 100*CoveredSize/SizeT1
            PercCover2 = 100*CoveredSize/SizeT2
            if PercCover1>Coverage or PercCover2>Coverage:
                Cover = True
    elif Beg>=Beg2 and End<=End2:
        Cover = True
    else:
        Cover = False
    return Cover
        
        
def split_per_coverage(GroupsMembers, DicoPos, Coverage):
    # This function distribute transcripts in groups according to their percentage of overlap
    FinalListe = [] 
    for Transcript in GroupsMembers:
        Data = DicoPos[Transcript]
        Beg = Data[1]
        End = Data[2]
        if len(FinalListe) == 0:
            FinalListe.append([Transcript])
        else:
            for OtherTransc in GroupsMembers:
                if OtherTransc!=Transcript:
                    Data2 = DicoPos[OtherTransc]
                    Beg2 = Data2[1]
                    End2 = Data2[2]
                    Cover = calculate_coverage(Beg, End, Beg2, End2, Coverage) ###
                    if Cover == False:
                        continue
                    else:
                        Transcripts1Present = False
                        Transcript2Present = False
                        AllCover = True
                        c = 0
                        PlacementDonne = False
                        for Liste in FinalListe:
                            PlaceDansListe = True
                            for Transcript3 in Liste:
                                if Transcript3 == Transcript:
                                    Transcripts1Present = True
                                elif Transcript3 == OtherTransc:
                                    Transcript2Present = True
                                else:
                                    Data3 = DicoPos[Transcript3]
                                    Beg3 = Data3[1]
                                    End3 = Data3[2]
                                    NewCover1 = calculate_coverage(Beg, End, Beg3, End3, Coverage)
                                    if NewCover1 == False:
                                        PlaceDansListe = False
                                        break
                                    NewCover2 = calculate_coverage(Beg2, End2, Beg3, End3, Coverage)
                                    if NewCover2 == False:
                                        PlaceDansListe = False
                                        break
                            if PlaceDansListe == True:
                                PlacementDonne = True
                                if Transcripts1Present == False:
                                    FinalListe[c].append(Transcript)
                                if Transcript2Present == False:
                                    FinalListe[c].append(OtherTransc)
                            c+=1
                        if PlacementDonne == False:
                            FinalListe.append([Transcript, OtherTransc])
    return FinalListe

                                
def sort_per_coverage(DicoFinal1, DicoPos, Coverage):
    # This function sort transcripts according to the percentage of overlap
    DicoFinal2 = {}
    c = 1
    for i in DicoFinal1.keys():
        GroupsMembers = DicoFinal1[i]
        if len(DicoFinal1[i]) == 1:
            NewName = "Orthogroups"+str(c)
            DicoFinal2[NewName] = GroupsMembers
            c+=1
        else:
            NewListe = split_per_coverage(GroupsMembers, DicoPos, Coverage)
            for Listes in NewListe:
                NewName = "Orthogroups"+str(c)
                DicoFinal2[NewName] = Listes
                c+=1
    print ("There is a final number of "+str(len(DicoFinal2)))
    return DicoFinal2


def build_filal_file(D, NameFile):
    # This function builds the 2 finals orthogroups files
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
    Groups = open_file("Allgroups.rst")
    DicoGroups = store_groups(Groups)
    Info = open_file("AllDeNovo_Info")
    assess_transcripts_size(Info)
    DicoPos = build_dic_pos(Info)
    DicoFinal1 = rebuild_groups(DicoGroups, DicoPos)
    DicoFinal_SortedPyStart = sort_by_starting_point(DicoFinal1, DicoPos, 500)
    DicoFinal_SortedPyCoverage = sort_per_coverage(DicoFinal1, DicoPos, 70)
    build_filal_file(DicoFinal1, "DeNovoTranscriptsGroups_Rst_Brut")
    build_filal_file(DicoFinal_SortedPyCoverage, "DeNovoTranscriptsGroups_Rst1")
    build_filal_file(DicoFinal_SortedPyStart, "DeNovoTranscriptsGroups_Rst2")

