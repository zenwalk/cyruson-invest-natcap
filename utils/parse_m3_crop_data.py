import os, sys, time, numpy, glob
from osgeo import gdal
from osgeo import osr
from osgeo.gdalconst import *

def die(message):
    print message
    sys.exit(1)

def countValidNonZeroData(band):
    noDataValue = band.GetNoDataValue()
    array = band.ReadAsArray(0,0,band.XSize,band.YSize)
    nonzero = numpy.where(array != noDataValue)[0].size
    total = band.XSize*band.YSize
    return (nonzero,total,nonzero/float(total))



#for timing runs
startTime = time.time()
gdal.AllRegister()

PATH = '../tmp_data/175crops/'
for filename in glob.glob(os.path.join(PATH, '*.nc')):
    
    #Register all drivers at once
    
    print filename
    
    #open file
    dataset = gdal.Open(filename, GA_ReadOnly);
    
    if dataset is None:
        die('Could not open ' + filename)
    
    geotransform = dataset.GetGeoTransform()
    originX = geotransform[0]
    originY = geotransform[3]
    pixelWidth = geotransform[1]
    pixelHeight = geotransform[5]
    
    yieldBand = dataset.GetRasterBand(1)
    areaBand = dataset.GetRasterBand(2)
    
    
    nonzero, total, percent = countValidNonZeroData(yieldBand)
    
    #for y in range(band.YSize):
    #    scanline = band.ReadAsArray(0,y,band.XSize, 1)
    #    nonzero += numpy.where(scanline != noDataValue)[0].size
    #    total += band.XSize
    
    #print 'RasterCount: ' + str(dataset.RasterCount)
    #print 'upper left: (' + str(originX) + ', ' + str(originY) + ')'
    #print 'pixel dimensions (' + str(pixelWidth) + ', ' + str(pixelHeight) + ')'
    
    print 'nonzero elements: ' + str(nonzero) + ' total: ' + str(total) + ' percent: ' + str(float(nonzero)/total)
    print 'Time elapsed: ' + str(time.time()-startTime) + ' seconds'
    print