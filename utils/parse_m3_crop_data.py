import os, sys, time, numpy, glob, re
from osgeo import gdal
from osgeo import osr
from osgeo.gdalconst import *

BANDIDS = {'yield':0, 'harvestArea':1}
MAXBANDS = len(BANDIDS)
PATH = '../tmp_data/175crops/'

def createCropIdMap(cropNames):
    """Build up cropId table"""
    cropIds = {}
    for cropName in cropNames: 
        if cropName not in cropIds:
            cropIds[cropName] = len(cropIds)
    return cropIds

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

#for timing runs
startTime = time.time()

#Register all drivers at once
gdal.AllRegister()

globalMap = {}
filenames = glob.glob(os.path.join(PATH, '*.nc'))
cropNameRe = re.compile("^(.*/)([^/]*)_globalMap = {}5min")
cropNames = [cropNameRe.match(filename).group(2) for filename in filenames]
cropIds = createCropIdMap(filenames)

for filename in filenames:
    
    #break the cropname (apple) out of the path /home/joe/path/apple_5min.nc
    cropName = cropNameRe.match(filename).group(2) 
    
    #create a cropId if that crop hasn't been seen before
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
    
    addValidCropDataToMap(yieldBand, cropIds[cropName], globalMap, BANDIDS['yield'])
    addValidCropDataToMap(harvestedAreaBand, cropIds[cropName], globalMap, BANDIDS['harvestArea'])
    
    print
    sys.stdout.flush()
