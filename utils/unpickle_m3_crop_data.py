#The purpose of this script is to unpickle the sparse M3 crop data and convert
#into a single coherent data structure that allows indexes by coordinate
#to map to a map of cropids which in turn map to average yields and fractional
#land usage.

import os, sys, time, numpy, glob, re, struct
from osgeo import gdal
from osgeo import osr
from osgeo.gdalconst import *

FILE = './globalMap_bin.pkl'

#Functions

def die(message):
    """Error handler"""
    print message
    sys.exit(1)

#main script

#This will map global coordinate tuples to a map of cropIds to yield/area arrays
globalMap = {}


def unpickleBinaryMap(filename):
    with open(filename,'rb') as inFile:
        nCoords = struct.unpack('!i',inFile.read(4))[0] # ncoords
        print nCoords
        globalMap = {}
        for i in range(nCoords):
            coord = (struct.unpack_from('!h',inFile.read(2))[0],struct.unpack_from('!h',inFile.read(2))[0]) # x, y coord
            #print coord
            nCrops = struct.unpack_from('!b',inFile.read(1))[0] #n crops
            #print nCrops
            cropMap = {}
            for j in range(nCrops):
                cropId = struct.unpack_from('!b',inFile.read(1))[0] #cropId
                #print cropId 
                cropArray = [None,None]
                for k in range(2):
                    cropVal = struct.unpack_from('!f',inFile.read(4))[0]  # crop value (either area or yield) 
                    if cropVal != -1.0:
                        cropArray[k] = cropVal
                cropMap[cropId] = cropArray
                #print cropArray
            globalMap[coord] = cropMap
        #print globalMap
        return globalMap

globalMap = unpickleBinaryMap(FILE)
for coord in globalMap:
    for crop in  globalMap[coord]:
        print str(coord) + ', ' + str(crop) + ',' + str(globalMap[coord][crop])
sys.stdout.flush()
