#The purpose of this script is to unpickle the sparse M3 crop data and convert
#into a single coherent data structure that allows indexes by coordinate
#to map to a map of cropids which in turn map to average yields and fractional
#land usage.

import os, sys, time, numpy, glob, re, pickle
from osgeo import gdal
from osgeo import osr
from osgeo.gdalconst import *

FILE = './globalMap.pkl'

#Functions

def die(message):
    """Error handler"""
    print message
    sys.exit(1)

#main script

#This will map global coordinate tuples to a map of cropIds to yield/area arrays
globalMap = {}

inFile = open(FILE,'rb')
globalMap = pickle.load(inFile)
for coord in globalMap:
    for crop in  globalMap[coord]:
        print str(coord) + ', ' + str(crop) + ',' + str(globalMap[coord][crop])
sys.stdout.flush()
