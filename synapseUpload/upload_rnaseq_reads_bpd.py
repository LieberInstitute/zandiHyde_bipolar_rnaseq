# module load python/2.7.9
import os
import sys
import synapseclient;
from synapseclient import Wiki, File, Project, Folder

## login
syn = synapseclient.login()

## read in files to upload
dat = open('fastqFiles_LIBD_bipolar.txt', 'rU')
synFolder = "syn12207389"
for line in dat:
    line = line.strip()
    columns = line.split()
    leftRead = columns[4]
    rightRead = columns[5]
    fLeft = File(leftRead, parent = synFolder)
    syn.store(fLeft)
    fRight = File(rightRead, parent = synFolder)
    syn.store(fRight)