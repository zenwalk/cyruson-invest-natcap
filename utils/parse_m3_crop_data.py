import os, sys, time, numpy, glob, re
from osgeo import gdal
from osgeo import osr
from osgeo.gdalconst import *


bandIds = {'yield':0, 'harvestArea':1}
MAXBANDS = len(bandIds)
PATH = '../tmp_data/175crops/'

cropNameRe = re.compile("^(.*/)([^/]*)_5min")

globalMap = {}
cropIds = {}

#Build up cropId table
for filename in glob.glob(os.path.join(PATH, '*.nc')):
    cropName = cropNameRe.match(filename).group(2)
    if cropName not in cropIds:
        cropIds[cropName] = len(cropIds)
    print cropName

MAXCROPS = len(cropIds)

def die(message):
    print message
    sys.exit(1)

def countValidNonZeroData(band):
    
    noDataValue = band.GetNoDataValue()
    startTime = time.time()
    array = band.ReadAsArray(0, 0, band.XSize, band.YSize)
    print 'Read as array time elapsed: ' + str(time.time() - startTime) + ' seconds'
    startTime = time.time()
    nonzero = numpy.where((array != noDataValue) & (array != 0))[0].size
    print 'Count nonzero time elapsed: ' + str(time.time() - startTime) + ' seconds'
    total = band.XSize * band.YSize
    
    return (nonzero, total, nonzero / float(total))

def collectValidNonZeroData(band, cropId, map, bandId):
    noDataValue = band.GetNoDataValue()
    startTime = time.time()
    array2d = band.ReadAsArray(0, 0, band.XSize, band.YSize)
    
    nonzeroValidIndices = numpy.where((array2d != noDataValue) & (array2d != 0))
    values = array2d[nonzeroValidIndices]
    #nonzeroValidCoordinates = zip(nonzeroValidIndices[0],nonzeroValidIndices[1])
    #for coord in nonzeroValidCoordinates:
    for i in range(len(values)):
        coord = (nonzeroValidIndices[0][i], nonzeroValidIndices[1][i])
        if coord not in map:
            map[coord] = {} #default map of crops
        if cropId not in map[coord]:
            map[coord][cropId] = [None]*MAXBANDS #default array of band types for that crop at that coordinate
        map[coord][cropId][bandId] = values[i]
        
    print 'Collect ' + str(len(values)) + ' nonzero data for ' + str(bandId) + ' elapsed: ' + str(time.time() - startTime) + ' seconds'
    print 'Map length: ' + str(len(map))
#for timing runs
startTime = time.time()

#Register all drivers at once
gdal.AllRegister()




for filename in glob.glob(os.path.join(PATH, '*.nc')):
    
    cropName = cropNameRe.match(filename).group(2)
    if cropName not in cropIds:
        cropIds[cropName] = len(cropIds)
    print cropName
    
    #open file
    dataset = gdal.Open(filename, GA_ReadOnly);
    
    if dataset is None:
        die('Could not open ' + filename)
    
    geotransform = dataset.GetGeoTransform()
    originX = geotransform[0]
    originY = geotransform[3]
    pixelWidth = geotransform[1]
    pixelHeight = geotransform[5]
    
    harvestedAreaBand = dataset.GetRasterBand(1)
    yieldBand = dataset.GetRasterBand(2)
    
    collectValidNonZeroData(yieldBand, cropIds[cropName], globalMap, bandIds['yield'])
    collectValidNonZeroData(harvestedAreaBand, cropIds[cropName], globalMap, bandIds['harvestArea'])
    
    #nonzero, total, percent = countValidNonZeroData(yieldBand)
    #print 'yield nonzero elements: ' + str(nonzero) + ' total: ' + str(total) + ' percent: ' + str(float(nonzero)/total)

    #nonzero, total, percent = countValidNonZeroData(harvestedAreaBand)
    #print 'area nonzero elements: ' + str(nonzero) + ' total: ' + str(total) + ' percent: ' + str(float(nonzero)/total)
    
    print
    sys.stdout.flush()
