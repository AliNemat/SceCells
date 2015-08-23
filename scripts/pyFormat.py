### This scrip corrects formatting issue that NSight cannot solve

import glob, os, sys
import fileinput
import re
from shutil import copyfile, move

fileType = "*.h"
folder = "../src/srcGPU"

def lineReplace(line):
    # For unknown reasons, Nsight has formatting issues for lines that
    # contain operator() and __host__ . This function corrects that behavior.
    lineNew = re.sub(r'__host__\s+__device__','__host__ __device__', line)
    if( line != lineNew):
        print (lineNew)
    return lineNew

def replaceAll(file):
    print(file)
    fileNew = file + ".new"
    fileBackup = file + ".bak"
    copyfile(file, fileBackup)
    with open(file) as infile, open(fileNew, 'w') as outfile:
        for line in infile:
            line = lineReplace(line)
            outfile.write(line)
    os.remove(file)
    move(fileNew, file)

os.chdir(folder)
for file in glob.glob(fileType):
    replaceAll(file) 
    
