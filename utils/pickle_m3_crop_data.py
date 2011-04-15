#The purpose of this script is to read all the M3 crop data and convert
#into a single coherent data structure that allows indexes by coordinate
#to map to a map of cropids which in turn map to average yields and fractional
#land usage.

import os, sys, time, numpy, glob, re, struct, mmap
from osgeo import gdal
from osgeo import osr
from osgeo.gdalconst import *

#Global constants

BANDIDS = {'harvestArea':0, 'yield':1}
MAXBANDS = len(BANDIDS)
PATH = '../tmp_data/175crops/'
OUTFILE = '../tmp_data/m3pickle.bin'

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
        #nonzeroValidIndices are in the form row/col (lat/lng)
        coord = (nonzeroValidIndices[0][i], nonzeroValidIndices[1][i])
        if coord not in map:
            map[coord] = {} #default map of crops
        if cropId not in map[coord]:
            map[coord][cropId] = [None] * MAXBANDS #default array of band types for that crop at that coordinate
        map[coord][cropId][bandId] = values[i]

def verifyGeotransformData(filenames):
    originX = None
    originY = None
    pixelWidth = None
    pixelHeight = None
    xnorthup = None
    ynorthup = None
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
            xnorthup = geotransform[2]
            ynorthup = geotransform[4]
            
        #check main variables, they should all be the same for each dataset    
        if originX != geotransform[0]:
            die("originX not the same: " + str(originX) + ' ' + str(geotransform[0]))
        if originY != geotransform[3]:
            die("originY not the same: " + str(originY) + ' ' + str(geotransform[3]))
        if pixelWidth != geotransform[1]:
            die("pixelWidth not the same: " + str(pixelWidth) + ' ' + str(geotransform[1]))
        if pixelHeight != geotransform[5]:
            die("pixelHeight not the same: " + str(pixelHeight) + ' ' + str(geotransform[5]))
            
        print 'originX', originX
        print 'originY', originY
        print 'pixelWidth', pixelWidth
        print 'pixelHeight', pixelHeight
        print 'xnorthup', xnorthup
        print 'ynorthup', ynorthup

def pickleBinaryMap(cropMap, geoTransform, cropIds, fileName):
    def write(format, var, file):
        """Simplifies writing to a file using struct compression"""
        file.write(struct.pack(format, var))

    with open(fileName, 'wb') as file:
        #file = mmap.mmap(f.fileno(),0)

        #dump the list of crops and their ids
        write('!b', len(cropIds), file) # number of crop ids
        for cropName in cropIds:
            write('!b',len(cropName),file) #length of the crop name
            write('!'+str(len(cropName))+'s',cropName,file) #the crop name
            write('!b',cropIds[cropName],file)
        
        #write geotransform
        #write lngOrig, latOrig, lngWidth, latWidth  - from http://www.gdal.org/gdal_tutorial.html
        for v in [geoTransform[0], geoTransform[3], geoTransform[1], geoTransform[5]]:
            write('!f', v, file)

        write('!i', len(cropMap), file) # ncoords
        for coord in cropMap: 
            write('!h', coord[0], file) #latIndex
            write('!h', coord[1], file) #lngIndex
            write('!B', len(cropMap[coord]), file) # n crops 
            for cropId in cropMap[coord]:
                write('!B', cropId, file) # cropId
                for val in cropMap[coord][cropId]: # crop value (either area or yield)
                    write('!f', val if val != None else - 1.0, file)

#main script

#Register all drivers at once
gdal.AllRegister()

#This will map global coordinate tuples to a map of cropIds to yield/area arrays
globalMap = {}

#regular expression to break the cropname out of the path 
#example (apple) from /home/joe/path/apple_5min.nc
cropNameRe = re.compile("^(.*/)([^/]*)_5min")
cropIds = {}


filenames = glob.glob(os.path.join(PATH, 'banana*.nc'))

#verify that all geotransform data is the same across all files
print 'Verify geotransform data...',
verifyGeotransformData(filenames)
print 'passed!'

totalRuns = len(filenames)
runNumber = 0
totalTime = 0

geoTransform = None

for filename in filenames:
    
    if runNumber >= totalRuns:
        break
    runNumber += 1
    
    #Group 1 has the path, group 2 will have the filename in it 
    cropName = cropNameRe.match(filename).group(2) 
    
    #create a cropId if that crop hasn't been seen before
    if cropName not in cropIds:
        cropIds[cropName] = len(cropIds)
        #if cropName == 'banana':
        #    print cropName, cropIds[cropName]
        #    sys.exit(0)
    
    #continue
    
    print cropName + ': ' + str(runNumber) + ' of ' + str(totalRuns),

    #for timing runNumber
    startTime = time.time()
    
    #open file
    dataset = gdal.Open(filename, GA_ReadOnly);
    
    if dataset is None:
        die('Could not open ' + filename)
    
    #Get geotransform data which will later be saved
    geoTransform = dataset.GetGeoTransform()
    
    harvestedAreaBand = dataset.GetRasterBand(1)
    yieldBand = dataset.GetRasterBand(2)
    
    addValidCropDataToMap(yieldBand, cropIds[cropName], globalMap, BANDIDS['yield'])
    addValidCropDataToMap(harvestedAreaBand, cropIds[cropName], globalMap, BANDIDS['harvestArea'])
    
    del dataset # cause python garbage collection to deallocate
    
    currentTime = time.time() - startTime
    totalTime += currentTime
    
    print ' time: ' + str(currentTime) 
    sys.stdout.flush()
print "read time , " + str(totalTime)
pickleBinaryMap(globalMap, geoTransform, cropIds, OUTFILE)
print "pickle time , " + str(time.time()-startTime)
sys.stdout.flush()
