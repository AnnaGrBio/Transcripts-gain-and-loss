import os
os.chdir("path to folder")
import random


""" docstring

    This program builds the first info file with information about transcripts (splicing, size, etc).
    
    input
    ------------
    - gtf file of transcriptome assembly
    - renamed assembled transcripts
    
    output
    ------------
    - information file
    
"""


def open_file(NameFile):
    F = open(NameFile, "r")
    L = F.readlines()
    return L  


def nb_exon_per_transcript(File):
    # This function extracts the number of exons in transcripts from the gtf file of the transcriptome assembly
    # It returns a dictionary
    Dico = {}
    transcript_name = ""
    nb_exon = 0
    for i in File[2:]:
        ligne = i.split()
        if ligne[2] == "transcript":
            if transcript_name!="":
                Dico[transcript_name] = str(nb_exon)
                nb_exon = 0
            transcript_name = ligne[11][1:len(ligne[11])-2]
        else:
            nb_exon+=1
    Dico[transcript_name] = str(nb_exon)
    return Dico


def calculate_spliced_size(File):
    # This function calculates the size of the spliced transcripts
    Dico={}
    for i in File:
        if i[0]==">":
            #print i
            transcript_name = i.split("_")[0][1:]
            #print transcript_name
        else:
            ligne = i.split("\n")[0]
            Size = len(ligne)
            Dico[transcript_name] = str(Size)
    return Dico


def fill_dic_combined_information(File, DicoExons, DicoSplicedSize, FinalFileName):
    # This function combine the two dictionary of transcript size and nb of exons with all other informations contained in the gtf file, 
    # and creates a final output information file.
    D = {}
    for i in File[2:]:
        ligne = i.split()
        if ligne[2] == "transcript":
            Chrom = ligne[0]
            position_begin = ligne[3]
            position_end = ligne[4]
            direction = ligne[6]
            gene_name =ligne[9][1:len(ligne[9])-2]
            transcript_name = ligne[11][1:len(ligne[11])-2]
            Totaltranscript_name = transcript_name+"_"+gene_name
            Cov = ligne[13][1:len(ligne[13])-2]
            FPKM = ligne[15][1:len(ligne[15])-2]
            TPM = ligne[17][1:len(ligne[17])-2]
            UnsplicedSize = int(position_end)-int(position_begin)+1
            UnsplicedSize = str(UnsplicedSize)
            SplicedSize = DicoSplicedSize[transcript_name]
            nb_exon = DicoExons[transcript_name]
            Liste = [Chrom, position_begin, position_end, direction, Totaltranscript_name, transcript_name, UnsplicedSize, SplicedSize, nb_exon, Cov, FPKM, TPM]
            if gene_name not in D:
                D[gene_name] = [Liste]
            else:
                D[gene_name].append(Liste)
    F = open(FinalFileName, "w")
    F.write("gene_name"+","+"NbTranscripts"+","+"Chrom"+","+"position_begin"+","+"position_end"+","+"direction"+","+"Totaltranscript_name"+","+"transcript_name"+","+"UnsplicedSize"+","+"SplicedSize"+","+"nb_exon"+","+"Cov"+","+"FPKM"+","+"TPM"+"\n")
    for i in D.keys():
        NbTranscripts = str(len(D[i]))
        NameGene = i
        for j in D[i]:
            F.write(NameGene+","+NbTranscripts)
            for k in j:
                F.write(","+k)
            F.write("\n")
    F.close()
            
        
def main_function():
    # Has to be reshuffled for each line
    File = open_file("AK5stringtie.gtf")
    DicoExons = nb_exon_per_transcript(File)
    FileSpliced = open_file("AK5RenamedTranscripts")
    DicoSplicedSize = calculate_spliced_size(FileSpliced)
    fill_dic_combined_information(File, DicoExons, DicoSplicedSize, "AK5_InfoFile")




