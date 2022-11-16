import cv2
import os
os.chdir("path to folder")
import random


""" docstring

    This program attribute a genomic position to each de novo transcript
    
    input
    ------------
    - bed file of de novo transcripts
    
    output
    ------------
    - position in genome of de novo transcripts
"""


def open_file(NameFile):
    F=open(NameFile, "r")
    L=F.readlines()
    return L  


def access_exonic_category(NewName, Signe):
    # access genomic position of de novo transcripts
    F = open(NewName, "r")
    L = F.readlines()
    CoordBegTransc = int(L[0].split()[1])
    CoordEndTranscr = int(L[0].split()[2])
    print CoordBegTransc
    print CoordEndTranscr
    ListeAllPosBegin = []
    ListeAllPosEnd = []
    for i in L[1:]:
        SigneOfOfficialGene = i.split()[10]
        if SigneOfOfficialGene == Signe:
            RefPosBeg = int(i.split()[6])
            RefPosEnd = int(i.split()[7])
            ListeAllPosBegin.append(RefPosBeg)
            ListeAllPosEnd.append(RefPosEnd)
    BeginEarlier = True
    for i in ListeAllPosBegin:
        if i<CoordBegTransc:
            BeginEarlier = False
            break
    EndLater = True
    for i in ListeAllPosEnd:
        if i>CoordEndTranscr:
            EndLater = False
            break
    if BeginEarlier == True or EndLater == True:
        return "ExonLonger"
    elif BeginEarlier == True and EndLater == True:
        return "ExonLonger"
    else:
        return "ExonInside"
        
        
def access_._category(NewName, Signe):
    # Get undirectional exon finaly forward
    F = open(NewName, "r")
    L = F.readlines()
    CoordBegTransc = int(L[0].split()[1])
    CoordEndTranscr = int(L[0].split()[2])
    ListeAllPosBegin = []
    ListeAllPosEnd = []
    for i in L[1:]:
        SigneOfOfficialGene = i.split()[10]
        RefPosBeg = int(i.split()[6])
        RefPosEnd = int(i.split()[7])
        ListeAllPosBegin.append(RefPosBeg)
        ListeAllPosEnd.append(RefPosEnd)
    BeginEarlier = True
    for i in ListeAllPosBegin:
        if i<CoordBegTransc:
            BeginEarlier = False
            break
    EndLater = True
    for i in ListeAllPosEnd:
        if i>CoordEndTranscr:
            EndLater = False
            break
    if BeginEarlier == True or EndLater == True:
        return "UndirectionalExonLonger"
        
    else:
        return "UndirectionalExonInside"


def define_category(NewName,Signe):
    # determine position
    if Signe == ".":
        Signe = "+"
    F = open(NewName, "r")
    L = F.readlines()
    FinalCategory = ""
    if len(L) == 1:
        FinalCategory = "Intergenic"
    else:
        Elements = []
        TheSignesAreTheSame = False
        for i in L[1:]:
            SigneOfOfficialGene = i.split()[10]
            if SigneOfOfficialGene == Signe:
                TheSignesAreTheSame = True
                TypeOfEltRecovered = i.split()[12]
                if TypeOfEltRecovered not in Elements:
                    Elements.append(TypeOfEltRecovered)
            else:
                TypeOfEltRecovered = i.split()[12]
                if TypeOfEltRecovered not in Elements:
                    Elements.append(TypeOfEltRecovered)
                
        if TheSignesAreTheSame == False and Signe!=".":
            if "ncRNA_gene" in Elements or "Pseudogene" in Elements and "CDS" not in Elements:
                FinalCategory = "Intergenic"
            else:
                FinalCategory = "ReverseGenic"
                print Elements
        elif TheSignesAreTheSame == False and Signe==".":
            FinalCategory = access_._category(NewName, Signe)
        else:
            if "pseudogene" in Elements and "CDS" not in Elements: 
                FinalCategory = "Pseudogene"
            elif "ncRNA_gene" in Elements and "CDS" not in Elements:
                FinalCategory = "NcRNA"
            elif "gene" in Elements and "exon" not in Elements:
                FinalCategory = "Intronic"
            elif "gene" in Elements and "exon" in Elements:
                FinalCategory = access_exonic_category(NewName, Signe)
            elif "gene" not in Elements and "exon" in Elements:
                FinalCategory = access_exonic_category(NewName, Signe)
            else:
                FinalCategory = "Weird"
    return FinalCategory


def run_all_functions(File, NameFinalFile):
    c = 0
    DicoFinal = {}
    for i in File:
        Signe = i.split("	")[3]
        NameORF = i.split("	")[4].split("\n")[0]
        Intermediar = open("Intermediar", "w")
        Intermediar.write(i)
        Intermediar.close()
        NewName = "Test"   #+str(c)
        os.system('bedtools intersect -wa -wb -a Intermediar -b DrosoGFF.bed -sorted > ' +NewName)
        Category = define_category(NewName, Signe)
        if Category == "ReverseGenic":
            print NameORF
        if Category == "Weird":
            break
        DicoFinal[NameORF] = Category
        
        c+=1
    F = open(NameFinalFile, "w")
    for i in DicoFinal.keys():
        F.write(i+","+DicoFinal[i]+"\n")
    F.close()


def run_all_functions_function():
    Data = open_file("AK5_bed")
    run_all_functions(Data, "AK5_PositionsInGenome")
