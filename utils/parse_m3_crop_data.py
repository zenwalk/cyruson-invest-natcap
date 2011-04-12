import os, sys, time, numpy, glob, re
from osgeo import gdal
from osgeo import osr
from osgeo.gdalconst import *

def die(message):
    print message
    sys.exit(1)

def countValidNonZeroData(band):
    
    noDataValue = band.GetNoDataValue()
    startTime = time.time()
    array = band.ReadAsArray(0,0,band.XSize,band.YSize)
    print 'Read as array time elapsed: ' + str(time.time()-startTime) + ' seconds'
    startTime = time.time()
    nonzero = numpy.where((array != noDataValue) & (array != 0))[0].size
    print 'Count nonzero time elapsed: ' + str(time.time()-startTime) + ' seconds'
    total = band.XSize*band.YSize
    
    return (nonzero,total,nonzero/float(total))

#for timing runs
startTime = time.time()

#Register all drivers at once
gdal.AllRegister()

PATH = '../tmp_data/175crops/'

cropNameRe = re.compile("^(.*/)([^/]*)_5min")

for filename in glob.glob(os.path.join(PATH, '*.nc')):
    
    cropName = cropNameRe.match(filename).group(2)
    print cropName
    
    #open file
    startTime = time.time()
    dataset = gdal.Open(filename, GA_ReadOnly);
    print 'Open file time elapsed: ' + str(time.time()-startTime) + ' seconds'
    
    if dataset is None:
        die('Could not open ' + filename)
    
    geotransform = dataset.GetGeoTransform()
    originX = geotransform[0]
    originY = geotransform[3]
    pixelWidth = geotransform[1]
    pixelHeight = geotransform[5]
    
    startTime = time.time()
    harvestedAreaBand = dataset.GetRasterBand(1)
    yieldBand = dataset.GetRasterBand(2)
    print 'Get bands time elapsed: ' + str(time.time()-startTime) + ' seconds'
    
    
    nonzero, total, percent = countValidNonZeroData(yieldBand)
    print 'yield nonzero elements: ' + str(nonzero) + ' total: ' + str(total) + ' percent: ' + str(float(nonzero)/total)

    nonzero, total, percent = countValidNonZeroData(harvestedAreaBand)
    print 'area nonzero elements: ' + str(nonzero) + ' total: ' + str(total) + ' percent: ' + str(float(nonzero)/total)
    
    print