#The purpose of this script is to unpickle the sparse M3 crop data and convert
#into a single coherent data structure that allows indexes by coordinate
#to map to a map of cropids which in turn map to average yields and fractional
#land usage.

import os, sys, time, numpy, glob, re, struct, mmap

FILENAME = './globalMap_bin.pkl'

#Functions

def die(message):
    """Error handler"""
    print message
    sys.exit(1)

def unpickleCoordinates(inFileMap):
    """Looks up coordinates and file positions for later access, this way
    we don't load the entire dataset at once, only as needed"""
    nCoords = struct.unpack('!i',inFileMap.read(4))[0] # ncoords
    print nCoords
    globalMap = {}
    for i in range(nCoords):
        if i % 10000 == 0: #follow progress
            percent = i/float(nCoords)
            print percent
            sys.stdout.flush()
        coord = (struct.unpack_from('!h',inFileMap.read(2))[0],struct.unpack_from('!h',inFileMap.read(2))[0]) # x, y coord
        nCrops = struct.unpack_from('!B',inFileMap.read(1))[0] #n crops
        #the -1 comes from the fact that we already read the nCrops, but we want to re-read it on the load
        globalMap[coord] = inFileMap.tell()-1
        #jump over the crop data 
        inFileMap.seek((1+4+4)*nCrops,os.SEEK_CUR)
    return globalMap

def getCropData(inFileMap,filePos):
    """Parses out the crop data for the given file position"""
    inFileMap.seek(filePos)
    nCrops = struct.unpack_from('!B',inFileMap.read(1))[0] #n crops
    
    cropMap = {}
    for j in range(nCrops):
        cropId = struct.unpack_from('!B',inFileMap.read(1))[0] #cropId
        print cropId 
        cropArray = [None,None]
        for k in range(2):
            cropVal = struct.unpack_from('!f',inFileMap.read(4))[0]  # crop value (either area or yield) 
            if cropVal != -1.0:
                cropArray[k] = cropVal
        cropMap[cropId] = cropArray
    return cropMap
    

#This will map global coordinate tuples to a map of file positions that can 
#later be read off disk to yield/area arrays
globalMap = {}
with open(FILENAME,'r+b') as inFile:
    inFileMap = mmap.mmap(inFile.fileno(),0)
    globalMap = unpickleCoordinates(inFileMap)

#globalMap is loaded, for fun look up some crop data    
for coord in globalMap:
    print coord, globalMap[coord], getCropData(inFileMap, globalMap[coord])
    keys = getCropData(inFileMap, globalMap[coord]).keys()
    keys.sort()
    print keys
    sys.stdout.flush()
    sys.exit(0)
