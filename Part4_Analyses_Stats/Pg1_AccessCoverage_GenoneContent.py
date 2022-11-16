import os
os.chdir("path to folder")
import random


""" docstring

    This program calculates the percentage of content of genomes elements
    
    input
    ------------
    - GTF file of D mel reference genome
    
    output
    ------------
    - Percentage of coverage of different elements in genomes
    
"""


def open_file(NameFile):
    F=open(NameFile, "r")
    L=F.readlines()
    return L  


#Data chromosomes sizes
DicoSizeChrom = {"2L":23513712, "2R":25286936, "3L":28110227, "3R":32079331, "4":1348131, "X":23542271, "Y":3667352, "mitochondrion_genome":19524}

# 1. Acess content in NcRNA and pseudogenes. The function uses the GFF file with all coordinates.
# When one of these elements in met, the coordinates are retreived, and the genome coverage is calculated.
def access_other_content(GFFfile): 
    TotalSum = (23513712+25286936+28110227+32079331+1348131+23542271+3667352+19524)*2
    TotalNcRNA = 0
    TotalPseudogene = 0
    c = 0
    for i in GFFfile:
        ligne = i.split("	")
        if len(ligne)>2:
            if ligne[1] == "FlyBase" and ligne[2] == "ncRNA_gene":
                SizeNcRNA = int(ligne[4])-int(ligne[3])
                TotalNcRNA += SizeNcRNA
                
            elif ligne[1] == "FlyBase" and ligne[2] == "pseudogene":
                SizePseudogene = int(ligne[4])-int(ligne[3])
                TotalPseudogene += SizePseudogene
            
    PercentageNcRNA = float(100)*float(TotalNcRNA)/float(TotalSum)
    print ("Percentage NcRNA content")
    print (PercentageNcRNA)
    print ("***********************") 
    PercentagePseudogene = float(100)*float(TotalPseudogene)/float(TotalSum)
    print ("Percentage Pseudogene content")
    print (PercentagePseudogene)
    print ("***********************")
    return PercentageNcRNA,PercentagePseudogene
     

# 2. Acess content in Genes. The function uses the GFF file with all coordinates.
## When one of these elements in met, the coordinates are retreived, and the genome coverage is calculated.
## All genes coordinates are firt stored in a dictionary, and added to this dictionary by calculating if they overlap to a gene already present.
## If the new gene to be added is inside an existing gene, it is not added. If it starts before and end after an existing gene, this gene is added and the other one removed. If they just overlap, a new gene is built based on the first start and the last stop. All this according to the frame
## The dico is then resorted, and the genes coverage calculated.


def reduce_dic(Dico):
    # This function to sort the final dictionary of genes. Some might still still overlapp after the first round of sorting. As it does not make any difference as the sorting 1 works well enough, it is not used for this dataset.
    for Chrom in Dico.keys():
        ListeGeneChrom = Dico[Chrom]
        for Genes in ListeGeneChrom:
            GeneStart = Genes[0]
            GeneEnd = Genes[1]
            for OtherGene in ListeGeneChrom:
                StartOtherGene = OtherGene[0]
                EndOtherGene = OtherGene[1]
                if OtherGene == Genes:
                    continue
                else:
                    if GeneStart<=StartOtherGene and GeneEnd>=EndOtherGene:
                            ListeGeneChrom.remove(OtherGene)
                            break
                    elif GeneStart>=StartOtherGene and GeneEnd<=EndOtherGene:
                            ListeGeneChrom.remove(Genes)
                            break
                    elif GeneStart<StartOtherGene and GeneEnd<=EndOtherGene and GeneEnd>=StartOtherGene:
                            ListeGeneChrom.remove(Genes)
                            ListeGeneChrom.remove(OtherGene)
                            ListeGeneChrom.append([GeneStart,EndOtherGene])
                            break
                    elif GeneStart>=StartOtherGene and GeneEnd>EndOtherGene and GeneStart<=EndOtherGene:
                            ListeGeneChrom.remove(Genes)
                            ListeGeneChrom.remove(OtherGene)
                            ListeGeneChrom.append([StartOtherGene,GeneEnd])
                            break
        Dico[Chrom] = ListeGeneChrom
    return Dico


def fill_dic(DicoGenes,Chrom,GeneStart,GeneEnd):
    # This function fill the dictionaries according to overlaps
    ListeGeneChrom = DicoGenes[Chrom]
    if len(ListeGeneChrom) == 0:
        DicoGenes[Chrom].append([GeneStart,GeneEnd])
    else:
        Action = False
        for Genes in ListeGeneChrom:
            StartOtherGene = Genes[0]
            EndOtherGene = Genes[1]
            if GeneStart<=StartOtherGene and GeneEnd>=EndOtherGene:
                ListeGeneChrom.remove(Genes)
                ListeGeneChrom.append([GeneStart,GeneEnd])
                DicoGenes[Chrom] = ListeGeneChrom
                Action = True
                break
            elif GeneStart>=StartOtherGene and GeneEnd<=EndOtherGene:
                Action = True
                break
            elif GeneStart<StartOtherGene and GeneEnd<=EndOtherGene and GeneEnd>=StartOtherGene:
                ListeGeneChrom.remove(Genes)
                ListeGeneChrom.append([GeneStart,EndOtherGene])
                DicoGenes[Chrom] = ListeGeneChrom
                Action = True
                break
            elif GeneStart>=StartOtherGene and GeneEnd>EndOtherGene and GeneStart<=EndOtherGene:
                ListeGeneChrom.remove(Genes)
                ListeGeneChrom.append([StartOtherGene,GeneEnd])
                DicoGenes[Chrom] = ListeGeneChrom
                Action = True
                break
        if Action == False:
            DicoGenes[Chrom].append([GeneStart,GeneEnd])
    return DicoGenes


def access_gene_content(GFFfile):
    ## Main function
    TotalSum = (23513712+25286936+28110227+32079331+1348131+23542271+3667352+19524)*2
    TotalGene = 0
    DicoGenesForward = {"2L":[], "2R":[], "3L":[], "3R":[], "4":[], "X":[], "Y":[], "mitochondrion_genome":[]}
    DicoGenesReverse = {"2L":[], "2R":[], "3L":[], "3R":[], "4":[], "X":[], "Y":[], "mitochondrion_genome":[]}
    DicoGenesNames = {}
    TotalNbGenes = 0
    for i in GFFfile:
        ligne = i.split("	")
        if len(ligne)>2:
            if ligne[1] == "FlyBase" and ligne[2] == "gene":
                GeneStart = int(ligne[3])
                GeneEnd = int(ligne[4])
                GeneName = (ligne[8].split(";")[0]).split(":")[1]
                DicoGenesNames[GeneName]=[GeneStart,GeneEnd]
                Chrom = ligne[0]
                Signe = ligne[6]
                if Signe == "+":
                    TotalNbGenes+=1
                    DicoGenesForward = fill_dic(DicoGenesForward,Chrom,GeneStart,GeneEnd)
                else:
                    TotalNbGenes+=1
                    DicoGenesReverse = fill_dic(DicoGenesReverse,Chrom,GeneStart,GeneEnd)
    DicoGeneTotal = {"2L":[], "2R":[], "3L":[], "3R":[], "4":[], "X":[], "Y":[], "mitochondrion_genome":[]}
    for i in DicoGenesForward.keys():
        TheList = DicoGenesForward[i]
        for j in TheList:
            DicoGeneTotal[i].append(j)
    for i in DicoGenesReverse.keys():
        Liste = DicoGenesReverse[i]
        for j in Liste:
            DicoGeneTotal[i].append(j)
    CompteurGenesForward = 0
    for i in DicoGenesForward.keys():
        CompteurGenesForward+=(len(DicoGenesForward[i]))
    CompteurGenesReverse = 0
    for i in DicoGenesReverse.keys():
        CompteurGenesReverse+=(len(DicoGenesReverse[i]))
    CompteurGenes = 0
    for i in DicoGeneTotal.keys():
        CompteurGenes+=(len(DicoGeneTotal[i]))
    for Chrom in DicoGeneTotal.keys():
        ListeGenes = DicoGeneTotal[Chrom]
        for gene in ListeGenes:
            TotalGene+= (gene[1]-gene[0])
    PercentageGene = float(100)*float(TotalGene)/float(TotalSum)
    print ("Percentage genes content")
    print (PercentageGene)
    print ("***********************")
    FinalListeGenes = []
    for Chrom in DicoGeneTotal.keys():
        ListeGenes = DicoGeneTotal[Chrom]
        for Coordgene in ListeGenes:
            for GeneNames in DicoGenesNames.keys():
                if DicoGenesNames[GeneNames] == Coordgene:
                    FinalListeGenes.append(GeneNames)
                    break
    return FinalListeGenes,PercentageGene ### Return the name of the genes that are selected, so without overlap
                

# 3. Access content in introns. This is the most complicated, because the coordinates of the introns has to be built. For each gene of our list of genes, I first search which transcript is the longest, and select it. Then, i access all exons presents in the transcripts. The cool thing is that they are already sorted from first to last one in genomic positions. The i calculate the lengh of the introns inside, and calculate the percentage of the genome it represents.
# The biais is that some introns from transcript of one gene that would not overlap are not considered. But on the other hand, we find very high intron coverage, meaning that they are longer than exons, which i did not expect i have to search more
def access_intron_content(GFFfile, ListeGenes):
    TotalSum = (23513712+25286936+28110227+32079331+1348131+23542271+3667352+19524)*2
    DicoGene_LongerTranscript = {}
    TotalIntron = 0
    for i in GFFfile:
        ligne = i.split("	")
        if len(ligne)>7:
            Data = ligne[8]
            GeneOrTranscript = Data.split(":")[0]
            if GeneOrTranscript == "ID=transcript" and ligne[2] == "mRNA":
                SizeTranscript = int(ligne[4])-int(ligne[3])
                NomTranscript = (ligne[8].split(";")[0]).split(":")[1]
                NomGene = (ligne[8].split(";")[1]).split(":")[1]
                if NomGene in ListeGenes: ### Ici shift
                    if NomGene not in DicoGene_LongerTranscript.keys():
                        DicoGene_LongerTranscript[NomGene] = [NomTranscript,SizeTranscript]
                    else:
                        if DicoGene_LongerTranscript[NomGene][1]<SizeTranscript:
                            DicoGene_LongerTranscript[NomGene] = [NomTranscript,SizeTranscript]
    DicoTranscripts = {}
    for i in DicoGene_LongerTranscript.keys():
        NomTranscript = DicoGene_LongerTranscript[i][0]
        DicoTranscripts[NomTranscript] = []
    for i in GFFfile:
        ligne = i.split("	")
        if len(ligne)>7 and ligne[2] == "exon":
            NomParentTranscript = (ligne[8].split(";")[0]).split(":")[1]
            if NomParentTranscript in DicoTranscripts.keys():
                DicoTranscripts[NomParentTranscript].append([int(ligne[3]),int((ligne[4]))])
    for Transcript in DicoTranscripts.keys():
        ListeExons = DicoTranscripts[Transcript]
        if len(ListeExons)>1:
            if len(ListeExons) == 2:
                TotalIntron += ListeExons[1][0]-ListeExons[0][1]
            for i in range(1,len(ListeExons)):
                TheExon = ListeExons[i]
                PreviousExon = ListeExons[i-1]
                TotalIntron += TheExon[0]-PreviousExon[1]
    PercentageIntron = float(100)*float(TotalIntron)/float(TotalSum)
    print ("Percentage introns content")
    print (PercentageIntron)
    print ("***********************")
    return PercentageIntron
                

def main_function():
    GFFfile = open_file("Drosophila_melanogaster.BDGP6.32.54.chr.gff3")
    FinalListeGenes,PercentageGene = access_gene_content(GFFfile)
    PercentageNcRNA,PercentagePseudogene = access_other_content(GFFfile)
    PercentageIntron = access_intron_content(GFFfile,FinalListeGenes)
    PercentageNonCoding = float(100)-PercentageGene-PercentageNcRNA-PercentagePseudogene
    print ("Percentage intergenic region content")
    print (PercentageNonCoding)
    print ("***********************")
