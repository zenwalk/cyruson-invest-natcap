#The purpose of this script is to read all the M3 crop data and convert
#into a single coherent data structure that allows indexes by coordinate
#to map to a map of cropids which in turn map to average yields and fractional
#land usage.

import os, sys, time, numpy, glob, re, struct
from osgeo import gdal
from osgeo import osr
from osgeo.gdalconst import *

#Global constants

BANDIDS = {'yield':0, 'harvestArea':1}
MAXBANDS = len(BANDIDS)
PATH = '../tmp_data/175crops/'

#Functions

def die(message):
    """Error handler"""
    print message
    sys.exit(1)

def addValidCropDataToMap(band, cropId, map, bandId):
    """Given the band for a particular cropId, find the nonzero valid
    data in that band and store it in the map under the particular bandId"""
    noDataValue = band.GetNoDataValue()
    array2d = band.ReadAsArray(0, 0, band.XSize, band.YSize)
    
    nonzeroValidIndices = numpy.where((array2d != noDataValue) & (array2d != 0))
    values = array2d[nonzeroValidIndices]
    for i in range(len(values)):
        coord = (nonzeroValidIndices[0][i], nonzeroValidIndices[1][i])
        if coord not in map:
            map[coord] = {} #default map of crops
        if cropId not in map[coord]:
            map[coord][cropId] = [None]*MAXBANDS #default array of band types for that crop at that coordinate
        map[coord][cropId][bandId] = values[i]

def verifyGeotransformData(filenames):
    originX = None
    originY = None
    pixelWidth = None
    pixelHeight = None
    nBands = None
    bandX = None
    bandY = None
    for filename in filenames:
        #open file
        dataset = gdal.Open(filename, GA_ReadOnly);
        
        if dataset is None:
            die('Could not open ' + filename)
        
        #Get geotransform data which will later be saved
        geotransform = dataset.GetGeoTransform()
        #initialize the main variables
        if originX == None:
            originX = geotransform[0]
            originY = geotransform[3]
            pixelWidth = geotransform[1]
            pixelHeight = geotransform[5]
            
        #check main variables, they should all be the same for each dataset    
        if originX != geotransform[0]:
            die("originX not the same: " + str(originX) + ' ' + str(geotransform[0]))
        if originY != geotransform[3]:
            die("originY not the same: " + str(originY) + ' ' + str(geotransform[3]))
        if pixelWidth != geotransform[1]:
            die("pixelWidth not the same: " + str(pixelWidth) + ' ' + str(geotransform[1]))
        if pixelHeight != geotransform[5]:
            die("pixelWidth not the same: " + str(pixelHeight) + ' ' + str(geotransform[5]))

def pickleBinaryMap(map,fileName):
    with open(fileName,'wb') as outFile:
        outFile.write(struct.pack('!i',len(map))) # ncoords
        for coord in map: 
            outFile.write(struct.pack('!h',coord[0])) # x coord
            outFile.write(struct.pack('!h',coord[1])) # y coord
            outFile.write(struct.pack('!b',len(map[coord]))) # n crops
            for cropId in map[coord]:
                outFile.write(struct.pack('!b',cropId)) # cropId
                for val in map[coord][cropId]: # crop value (either area or yield)
                    if val != None:
                        outFile.write(struct.pack('!f',val))
                    else:
                        outFile.write(struct.pack('!f',-1.0))


#main script



#Register all drivers at once
gdal.AllRegister()

#This will map global coordinate tuples to a map of cropIds to yield/area arrays
globalMap = {}

#regular expression to break the cropname out of the path 
#example (apple) from /home/joe/path/apple_5min.nc
cropNameRe = re.compile("^(.*/)([^/]*)_5min")
cropIds = {}


filenames = glob.glob(os.path.join(PATH, '*.nc'))

#verify that all geotransform data is the same across all files
print 'Verify geotransform data...',
verifyGeotransformData(filenames)
print 'passed!'

MAXRUNS = len(filenames)
runs= 0

totalTime = 0

for filename in filenames:
    
    if runs >= MAXRUNS:
        break
    runs += 1
    
    
    #Group 1 has the path, group 2 will have the filename in it 
    cropName = cropNameRe.match(filename).group(2) 
    
    #create a cropId if that crop hasn't been seen before
    if cropName not in cropIds:
        cropIds[cropName] = len(cropIds)
    
    print cropName + ': ' + str(runs) + ' of ' + str(MAXRUNS),

    #for timing runs
    startTime = time.time()
    
    #open file
    dataset = gdal.Open(filename, GA_ReadOnly);
    
    if dataset is None:
        die('Could not open ' + filename)
    
    #Get geotransform data which will later be saved
    geotransform = dataset.GetGeoTransform()
    originX = geotransform[0]
    originY = geotransform[3]
    pixelWidth = geotransform[1]
    pixelHeight = geotransform[5]
    
    harvestedAreaBand = dataset.GetRasterBand(1)
    yieldBand = dataset.GetRasterBand(2)
    
    addValidCropDataToMap(yieldBand, cropIds[cropName], globalMap, BANDIDS['yield'])
    addValidCropDataToMap(harvestedAreaBand, cropIds[cropName], globalMap, BANDIDS['harvestArea'])
    
    del dataset # cause python garbage collection to deallocate
    
    currentTime = time.time()-startTime
    totalTime += currentTime
    
    print ' time: '+ str(currentTime) 
    sys.stdout.flush()
print "totalTime , " + str(totalTime)
pickleBinaryMap(globalMap,'globalMap_bin.pkl')
sys.stdout.flush()
