#The purpose of this script is to unpickle the sparse M3 crop data and convert
#into a single coherent data structure that allows indexes by coordinate
#to map to a map of cropids which in turn map to average yields and fractional
#land usage.

import os, sys, time, numpy, glob, re, struct, mmap

class M3Loader:
    # This dictionary maps coordinates to file positions, or crop data if
    # those file positions have already been looked up
    #_globalMap = {}
    #_inFileMap = None

    def __init__(self,FILENAME):
        self._globalMap = {}
        with open(FILENAME,'r+b') as inFile:
            self._inFileMap = mmap.mmap(inFile.fileno(),0)
            self._unpickleCoordinates(self._inFileMap)
            print len(self._globalMap)
    
    def _unpickleCoordinates(self,inFileMap):
        """Looks up coordinates and file positions for later access, this way
        we don't load the entire dataset at once, only as needed"""
        nCoords = struct.unpack('!i',inFileMap.read(4))[0] # ncoords
        self._globalMap = {}
        for i in range(nCoords):
            if i % 100000 == 0: #follow progress
                percent = i/float(nCoords)
                print percent
                sys.stdout.flush()
            coord = (struct.unpack_from('!h',inFileMap.read(2))[0],struct.unpack_from('!h',inFileMap.read(2))[0]) # x, y coord
            nCrops = struct.unpack_from('!B',inFileMap.read(1))[0] #n crops
            #the -1 comes from the fact that we already read the nCrops, but we want to re-read it on the load
            self._globalMap[coord] = inFileMap.tell()-1
            #jump over the crop data 
            inFileMap.seek((1+4+4)*nCrops,os.SEEK_CUR)
    
    def _getCropData(self,filePos):
        """Parses out the crop data for the given file position"""
        self._inFileMap.seek(filePos)
        nCrops = struct.unpack_from('!B',self._inFileMap.read(1))[0] #n crops
        
        cropMap = {}
        for j in range(nCrops):
            cropId = struct.unpack_from('!B',self._inFileMap.read(1))[0] #cropId
            print cropId 
            cropArray = [None,None]
            for k in range(2):
                cropVal = struct.unpack_from('!f',self._inFileMap.read(4))[0]  # crop value (either area or yield) 
                if cropVal != -1.0:
                    cropArray[k] = cropVal
            cropMap[cropId] = cropArray
        return cropMap
    
    def getCropData(self,coord):
        """Report crop map index, if not cached, look it up"""
        keys = self._globalMap.keys()
        keys.sort()
        originX = -180.0
        originY = -90.0
        pixelWidth = 0.08333333 
        pixelHeight = 0.08333333
        for key in keys:
            print (originY+key[0]*pixelHeight,originX+key[1]*pixelWidth)
        #print len(self._globalMap)
        if self._globalMap[coord] == None:
            return None
        if type(self._globalMap[coord]) is int:
            self._globalMap[coord] = self._getCropData(self._globalMap[coord])
        return self._globalMap[coord]
        
if __name__ == '__main__':
    m3loader = M3Loader('../core_data/m3crop_invest.bin')        

    originX = -180.0
    originY = -90.0
    pixelWidth = 0.08333333 
    pixelHeight = 0.08333333
    
    x = -53.60
    y = -0.40
        
    j = int((x-originX)/pixelWidth)
    i = int((y-originY)/pixelHeight)
    
    #i=854
    #j=962
    print i,j
    print m3loader.getCropData((i,j))
    
        
    #globalMap is loaded, for fun look up some crop data    
#    for coord in globalMap:
#        print coord, globalMap[coord], getCropData(inFileMap, globalMap[coord])
#        cropData = getCropData(inFileMap, globalMap[coord])
#        for key in cropData:        
#            print cropData[key][0]
#            sys.stdout.flush()
#        sys.exit(0)
