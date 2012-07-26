from sets import Set
import os, sys
from numpy import *

class checkASCIIFile:
    def __init__(self, theFile):
        theFile = open(theFile, "r")
        self.line = theFile.readline()
        self.holder = []
        while self.line:
            self.holder.append(int(self.line.split(":")[0]))
            self.line = theFile.readline()
        theFile.close()

    def countDuplicatesInList(self, dupedList):
        uniqueSet = Set(item for item in dupedList)
        return [(item, dupedList.count(item)) for item in uniqueSet]
    
    def hasDuplicates(self):
        for x in self.countDuplicatesInList(self.holder): 
            if x[1] > 1:
                return True
        return False
    def isSorted(self):
        xhold = 0
        for x in self.holder:
            if xhold:
                if x < xhold:
                    return False
            xhold = x
        return True

def createfolders(theWorkspace):

    try:
        thefolders=["output","intermediate"]
        for folder in thefolders:
            if not os.path.exists(theWorkspace+"\\"+folder):
                os.mkdir(theWorkspace+"\\"+folder)
    except:
        print "Error creating folders"
        
def switchws(direction="up"):
    if direction=="down":
        gp.workspace=original_workspace+"\\intermediate"
    else:
        gp.workspace=original_workspace 


def readFileIntoArray(theFile, skip=0):
    matrixfile = open(theFile)
    line = matrixfile.readline()
    if skip:
        line = matrixfile.readline()
    distmax = []
    while line:
        myarr = line.split()
        distmax.append(myarr)
        line = matrixfile.readline()
    matrixfile.close()
    return array(distmax).astype(float)

#function counts duplicates in a list
#arguments: python list
#return: list with count of dups
def countDuplicatesInList(dupedList):
    uniqueSet = Set(item for item in dupedList)
    return [(item, dupedList.count(item)) for item in uniqueSet]

#function checks if list has duplicates
#arguments: python list from countDuplicatesInList
#return: boolean, true is has duplicates
def ListHasDuplicates(theList):
    for x in countDuplicatesInList(theList): 
        if x[1] > 1:
            return True
    return False
#function sorts a text file before its used for reclassing raster
def sortreclassfile(thefile):
    newlines = []
    filetosort = open(thefile)
    lines = filetosort.readlines()
    for line in lines:
        splat = line.split(":")
        splat = map(int, splat)
        newlines.append(splat)
    filetosort.close()
    filetosort = open(thefile, "w")
    newlines = sorted(newlines, key=lambda cover: cover[0])
    for line in newlines:
        filetosort.write(str(line[0])+":"+str(line[1])+"\n")
    filetosort.close()
