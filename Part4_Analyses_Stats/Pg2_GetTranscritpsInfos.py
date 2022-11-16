import os
os.chdir("path to folder")
import random


""" docstring

    This program format file that contain several informations concerning de novo transcripts: nb introns, nb splicing, TPM, cumulated size
    
    input
    ------------
    - All de novo info
    
    output
    ------------
    - For files formated with transcripts informations
    
"""


def open_file(NameFile):
    F = open(NameFile, "r")
    L = F.readlines()
    return L  


def get_average_dic(Dico):
    NewDico = {}
    for i in Dico.keys():
        Liste = Dico[i]
        Total = len(Liste)
        TheSum = sum(Liste)
        Average = float(TheSum)/float(Total)
        NewDico[i] = Average
    return NewDico


# Nb Introns per transcripts per regions
def make_final_intron_file(Dico1,Dico2,Dico3,Dico4,Dico5,Dico6,Dico7,NameFile):
    File = open(NameFile, "w")
    File.write("Pop,Position,NbIntron"+"\n")
    for Position in Dico1.keys():
        Liste = Dico1[Position]
        for Data in Liste:
            File.write("AK5,"+Position+","+str(Data)+"\n")
    for Position in Dico2.keys():
        Liste = Dico2[Position]
        for Data in Liste:
            File.write("DK5,"+Position+","+str(Data)+"\n")
    for Position in Dico3.keys():
        Liste = Dico3[Position]
        for Data in Liste:
            File.write("GI5,"+Position+","+str(Data)+"\n")
    for Position in Dico4.keys():
        Liste = Dico4[Position]
        for Data in Liste:
            File.write("SW5,"+Position+","+str(Data)+"\n")
    for Position in Dico5.keys():
        Liste = Dico5[Position]
        for Data in Liste:
            File.write("UM,"+Position+","+str(Data)+"\n")
    for Position in Dico6.keys():
        Liste = Dico6[Position]
        for Data in Liste:
            File.write("YE,"+Position+","+str(Data)+"\n")
    for Position in Dico7.keys():
        Liste = Dico7[Position]
        for Data in Liste:
            File.write("Zamb,"+Position+","+str(Data)+"\n")
    File.close()
            

def build_dic_region(File, NamePopSelected):
    DicoData = {}
    for i in File[1:]:
        ligne = i.split("\n")[0]
        ligne = ligne.split(",")
        NameGene = ligne[0]
        NamePop = NameGene.split("_")[0]
        if NamePop == NamePopSelected:
            NbIntron = int(ligne[10])-1
            Position = ligne[14]
            if Position not in DicoData.keys():
                DicoData[Position] = [NbIntron]
            else:
                DicoData[Position].append(NbIntron)
    DicoAverage = get_average_dic(DicoData)
    for i in DicoAverage.keys():
        print ("Category : "+str(i)+" = "+str(DicoAverage[i]))
        print ("***")
    return DicoData


def main_function_nb_introns():
    DataDeNovo = open_file("AllDeNovo_Info")
    print ("AK5")
    DicoAK5 = build_dic_region(DataDeNovo, "AK5")
    print ("***")
    print ("DK5")
    DicoDK5 = build_dic_region(DataDeNovo, "DK5")
    print ("***")
    print ("GI5")
    DicoGI5 = build_dic_region(DataDeNovo, "GI5")
    print ("***")
    print ("SW5")
    DicoSW5 = build_dic_region(DataDeNovo, "SW5")
    print ("***")
    print ("UM")
    DicoUM = build_dic_region(DataDeNovo, "UM")
    print ("***")
    print ("YE")
    DicoYE = build_dic_region(DataDeNovo, "YE")
    print ("***")
    print ("Zamb")
    DicoZamb = build_dic_region(DataDeNovo, "Zamb")
    print ("***")
    make_final_intron_file(DicoAK5,DicoDK5,DicoGI5,DicoSW5,DicoUM,DicoYE,DicoZamb,"IntronRsts.txt")
    

# Nb splicing per transcripts per regions
def make_final_file_splicing(Dico1,Dico2,Dico3,Dico4,Dico5,Dico6,Dico7,NameFile):
    File = open(NameFile, "w")
    File.write("Pop,Position,NbTranscripts"+"\n")
    for Position in Dico1.keys():
        Liste = Dico1[Position]
        for Data in Liste:
            File.write("AK5,"+Position+","+str(Data)+"\n")
    for Position in Dico2.keys():
        Liste = Dico2[Position]
        for Data in Liste:
            File.write("DK5,"+Position+","+str(Data)+"\n")
    for Position in Dico3.keys():
        Liste = Dico3[Position]
        for Data in Liste:
            File.write("GI5,"+Position+","+str(Data)+"\n")
    for Position in Dico4.keys():
        Liste = Dico4[Position]
        for Data in Liste:
            File.write("SW5,"+Position+","+str(Data)+"\n")
    for Position in Dico5.keys():
        Liste = Dico5[Position]
        for Data in Liste:
            File.write("UM,"+Position+","+str(Data)+"\n")
    for Position in Dico6.keys():
        Liste = Dico6[Position]
        for Data in Liste:
            File.write("YE,"+Position+","+str(Data)+"\n")
    for Position in Dico7.keys():
        Liste = Dico7[Position]
        for Data in Liste:
            File.write("Zamb,"+Position+","+str(Data)+"\n")
    File.close()


def get_max_position(Liste):
    DicoAllPos = {"NcRNA":0, "Intergenic":0, "ReverseGenic":0, "ExonLonger":0, "Intronic":0, "ExonInside":0, "Pseudogene":0}
    for i in Liste:
        DicoAllPos[i]+=1
    MaxPos = ""
    MaxValue = 0
    for Pos in DicoAllPos.keys():
        if DicoAllPos[Pos]>MaxValue:
            MaxValue = DicoAllPos[Pos]
            MaxPos = Pos
    return MaxPos


def sort_position(Dico):
    DicoFinal = {}
    for Gene in Dico.keys():
        ListePosition = Dico[Gene]
        if len(ListePosition) == 1:
            DicoFinal[Gene] = ListePosition[0]
        else:
            MainPosition = get_max_position(ListePosition)
            DicoFinal[Gene] = MainPosition
    return DicoFinal


def make_dic_gene_region(File,NamePopSelected):
    DicoNbTranscript = {}
    DicoPosition = {}
    for i in File[1:]:
        ligne = i.split("\n")[0]
        ligne = ligne.split(",")
        NameGene = ligne[0]
        NamePop = NameGene.split("_")[0]
        if NamePop == NamePopSelected:
            Position = ligne[14]
            if NameGene not in DicoNbTranscript.keys():
                DicoNbTranscript[NameGene] = 1
                DicoPosition[NameGene] = [Position]
            else:
                DicoNbTranscript[NameGene] += 1
                if Position not in DicoPosition[NameGene]:
                    DicoPosition[NameGene].append(Position)
    NewDicoPosition = sort_position(DicoPosition)
    FinalDico = {}
    for Gene in NewDicoPosition.keys():
        Position = NewDicoPosition[Gene]
        NbTranscripts = DicoNbTranscript[Gene]
        if Position not in FinalDico.keys():
            FinalDico[Position] = [NbTranscripts]
        else:
            FinalDico[Position].append(NbTranscripts)
    DicoAverage = get_average_dic(FinalDico)
    for i in DicoAverage.keys():
        print ("Category : "+str(i)+" = "+str(DicoAverage[i]))
        print ("***")
    return FinalDico


def get_nb_single_transcripts(Dico):
    Nb = 0
    for cat in Dico.keys():
        Nb+=len(Dico[cat])
    return Nb
    
                
def main_function_nb_splicing():
    DataDeNovo = open_file("AllDeNovo_Info")
    print ("AK5")
    DicoAK5 = make_dic_gene_region(DataDeNovo, "AK5")
    print ("***")
    print ("DK5")
    DicoDK5 = make_dic_gene_region(DataDeNovo, "DK5")
    print ("***")
    print ("GI5")
    DicoGI5 = make_dic_gene_region(DataDeNovo, "GI5")
    print ("***")
    print ("SW5")
    DicoSW5 = make_dic_gene_region(DataDeNovo, "SW5")
    print ("***")
    print ("UM")
    DicoUM = make_dic_gene_region(DataDeNovo, "UM")
    print ("***")
    print ("YE")
    DicoYE = make_dic_gene_region(DataDeNovo, "YE")
    print ("***")
    print ("Zamb")
    DicoZamb = make_dic_gene_region(DataDeNovo, "Zamb")
    print ("***")
    NbUniqTransc_AK5 = get_nb_single_transcripts(DicoAK5)
    NbUniqTransc_DK5 = get_nb_single_transcripts(DicoDK5)
    NbUniqTransc_GI5 = get_nb_single_transcripts(DicoGI5)
    NbUniqTransc_SW5 = get_nb_single_transcripts(DicoSW5)
    NbUniqTransc_UM = get_nb_single_transcripts(DicoUM)
    NbUniqTransc_YE = get_nb_single_transcripts(DicoYE)
    NbUniqTransc_Zamb = get_nb_single_transcripts(DicoZamb)
    print ("Total Nb Transcript uniq transcripts AK5 : "+str(NbUniqTransc_AK5))
    print ("Total Nb Transcript uniq transcripts DK5 : "+str(NbUniqTransc_DK5))
    print ("Total Nb Transcript uniq transcripts GI5 : "+str(NbUniqTransc_GI5))
    print ("Total Nb Transcript uniq transcripts SW5 : "+str(NbUniqTransc_SW5))
    print ("Total Nb Transcript uniq transcripts UM : "+str(NbUniqTransc_UM))
    print ("Total Nb Transcript uniq transcripts YE : "+str(NbUniqTransc_YE))
    print ("Total Nb Transcript uniq transcripts Zamb : "+str(NbUniqTransc_Zamb))
    make_final_file_splicing(DicoAK5,DicoDK5,DicoGI5,DicoSW5,DicoUM,DicoYE,DicoZamb,"SplicingRsts.txt")
   

# TPM per de novo transcripts vs established transcripts
def build_final_file_TPM(Dico1,Dico2,Dico3,Dico4,Dico5,Dico6,Dico7,NameFile):
    File = open(NameFile, "w")
    File.write("Pop,Position,TPM"+"\n")
    for Position in Dico1.keys():
        Liste = Dico1[Position]
        for Data in Liste:
            File.write("AK5,"+Position+","+str(Data)+"\n")
    for Position in Dico2.keys():
        Liste = Dico2[Position]
        for Data in Liste:
            File.write("DK5,"+Position+","+str(Data)+"\n")
    for Position in Dico3.keys():
        Liste = Dico3[Position]
        for Data in Liste:
            File.write("GI5,"+Position+","+str(Data)+"\n")
    for Position in Dico4.keys():
        Liste = Dico4[Position]
        for Data in Liste:
            File.write("SW5,"+Position+","+str(Data)+"\n")
    for Position in Dico5.keys():
        Liste = Dico5[Position]
        for Data in Liste:
            File.write("UM,"+Position+","+str(Data)+"\n")
    for Position in Dico6.keys():
        Liste = Dico6[Position]
        for Data in Liste:
            File.write("YE,"+Position+","+str(Data)+"\n")
    for Position in Dico7.keys():
        Liste = Dico7[Position]
        for Data in Liste:
            File.write("Zamb,"+Position+","+str(Data)+"\n")
    File.close()


def make_dico_TPM_region_de_novo(FileDeNovo,FileAllTranscripts,NamePopSelected):
    ListeNameGene = []
    DicoPosTPM = {"EstablishedTranscripts":[]}
    for i in FileDeNovo[1:]:
        ligne = i.split("\n")[0]
        ligne = ligne.split(",")
        NameGene = ligne[0]
        NamePop = NameGene.split("_")[0]
        if NamePop == NamePopSelected:
            Position = ligne[14]
            TPM = float(ligne[13])
            ListeNameGene.append(NameGene)
            if Position not in DicoPosTPM.keys():
                DicoPosTPM[Position] = [TPM]
            else:
                DicoPosTPM[Position].append(TPM)
    for i in FileDeNovo[1:]:
        ligne = i.split("\n")[0]
        ligne = ligne.split(",")
        NameGene = ligne[0]
        TPM = float(ligne[13])
        if NameGene not in ListeNameGene:
            DicoPosTPM["EstablishedTranscripts"].append(TPM)
    DicoAverage = get_average_dic(DicoPosTPM)
    for i in DicoAverage.keys():
        print ("Category : "+str(i)+" = "+str(DicoAverage[i]))
        print ("***")
    return DicoPosTPM
    
            
def main_function_TPM():
    DataDeNovo = open_file("AllDeNovo_Info")
    DataAllAK5 = open_file("AK5_New_InfoFile")
    DataAllDK5 = open_file("DK5_New_InfoFile")
    DataAllGI5 = open_file("GI5_New_InfoFile")
    DataAllSW5 = open_file("SW5_New_InfoFile")
    DataAllUM = open_file("UM_New_InfoFile")
    DataAllYE = open_file("YE_New_InfoFile")
    DataAllZamb = open_file("Zamb_New_InfoFile")
    print ("AK5")
    DicoAK5 = make_dico_TPM_region_de_novo(DataDeNovo,DataAllAK5, "AK5")
    print ("***")
    print ("DK5")
    DicoDK5 = make_dico_TPM_region_de_novo(DataDeNovo,DataAllDK5, "DK5")
    print ("***")
    print ("GI5")
    DicoGI5 = make_dico_TPM_region_de_novo(DataDeNovo,DataAllGI5, "GI5")
    print ("***")
    print ("SW5")
    DicoSW5 = make_dico_TPM_region_de_novo(DataDeNovo,DataAllSW5, "SW5")
    print ("***")
    print ("UM")
    DicoUM = make_dico_TPM_region_de_novo(DataDeNovo,DataAllUM, "UM")
    print ("***")
    print ("YE")
    DicoYE = make_dico_TPM_region_de_novo(DataDeNovo,DataAllYE, "YE")
    print ("***")
    print ("Zamb")
    DicoZamb = make_dico_TPM_region_de_novo(DataDeNovo,DataAllZamb, "Zamb")
    print ("***") 
    build_final_file_TPM(DicoAK5,DicoDK5,DicoGI5,DicoSW5,DicoUM,DicoYE,DicoZamb,"TPMRsts.txt")
            

# Cumulated Unspliced size
def access_cumulated_size(FileDeNovo,NamePopSelected):
    FinalSize = 0
    SizeTotalGenome = 137567484*2
    for i in FileDeNovo[1:]:
        ligne = i.split("\n")[0]
        ligne = ligne.split(",")
        NameGene = ligne[0]
        NamePop = NameGene.split("_")[0]
        if NamePop == NamePopSelected:
            UnsplicedSize = int(ligne[8])
            FinalSize+=UnsplicedSize
    Coverage = float(100)*float(FinalSize)/float(SizeTotalGenome)
    print (NamePopSelected+"Size : "+str(FinalSize)+"Percentage : "+str(Coverage))


def MainFunctionCumulatedSize():
    DataDeNovo = open_file("AllDeNovo_Info")
    access_cumulated_size(DataDeNovo,"AK5")
    access_cumulated_size(DataDeNovo,"DK5")
    access_cumulated_size(DataDeNovo,"GI5")
    access_cumulated_size(DataDeNovo,"SW5")
    access_cumulated_size(DataDeNovo,"UM")
    access_cumulated_size(DataDeNovo,"YE")
    access_cumulated_size(DataDeNovo,"Zamb")
    

def Main_function():        
    main_function_nb_introns()        
    main_function_nb_splicing()
    main_function_TPM()
    MainFunctionCumulatedSize()            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
