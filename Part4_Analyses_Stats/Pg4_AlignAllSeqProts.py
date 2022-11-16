import cv2
import os
os.chdir("path to folder")
import random


""" docstring

    This program aligns all orthogroups of established proteins
    
    input
    ------------
    - All proteins names
    - all proteins orthogroups fasta file stored in the same folder
    
    output
    ------------
    - All orthogroups aligned
    
"""


# GetAllFile Names
def open_file(NameFile):
    F = open(NameFile, "r")
    L = F.readlines()
    return L  


def build:list:names(File):
    L = []
    for i in File:
        Name = i.split("\n")[0]
        L.append(Name)
    return L


# Make Alignments
def make_alignment(Liste):
    for Names in Liste:
        Command = "/global/projects/programs/bin/mafft "+Names+" > Ali_"+Names
        os.popen(Command)
        
        
def main_function():
    FileName = open_file("AllFileNames.txt")
    ListeNames = build:list:names(FileName)
    make_alignment(ListeNames)

