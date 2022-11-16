import cv2
import os
os.chdir("path to folder")                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              
import random


""" docstring

        This program creates a new name for all de novo transcripts of all lines, in order to include the name of the line in the header
        
        input
        ------------
        - Info file of each lines
        
        output
        ------------
        - info file with updated names
""" 


def open_file(NameFile):
    F=open(NameFile, "r")
    L=F.readlines()
    return L  


def remane_with_line_name(File, Header):
    NewName = Header+"_New_InfoFile"
    F = open(NewName, "w")
    F.write(File[0])
    for i in File[1:]:
        ligne = i.split(",")
        L0 = Header+"_"+ ligne[0]
        L6 = Header+"_"+ ligne[6]
        L7 = Header+"_"+ligne[7]
        F.write(L0+","+ligne[1]+","+ligne[2]+","+ligne[3]+","+ligne[4]+","+ligne[5]+","+L6+","+L7+","+ligne[8]+","+ligne[9]+","+ligne[10]+","+ligne[11]+","+ligne[12]+","+ligne[13])


def main_function():
    File = open_file("AK5_InfoFile")
    remane_with_line_name(File, "AK5")
    File = open_file("DK5_InfoFile")
    remane_with_line_name(File, "DK5")
    File = open_file("GI5_InfoFile")
    remane_with_line_name(File, "GI5")
    File = open_file("SW5_InfoFile")
    remane_with_line_name(File, "SW5")
    File = open_file("UM_InfoFile")
    remane_with_line_name(File, "UM")
    File = open_file("YE_InfoFile")
    remane_with_line_name(File, "YE")
    File = open_file("Zamb_InfoFile")
    remane_with_line_name(File, "Zamb")
