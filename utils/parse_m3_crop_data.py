import os, sys, time, numpy
from osgeo import gdal
from osgeo import osr
from osgeo.gdalconst import *

def die(message):
    print message
    sys.exit(1)

#for timing runs
startTime = time.time()

FILENAME = '../tmp_data/175crops/abaca_5min.nc'

#Register all drivers at once
gdal.AllRegister()

#open file
dataset = gdal.Open(FILENAME, GA_ReadOnly);

if dataset is None:
    die('Could not open ' + FILENAME)

geotransform = dataset.GetGeoTransform()
originX = geotransform[0]
originY = geotransform[3]
pixelWidth = geotransform[1]
pixelHeight = geotransform[5]

yieldBand = dataset.GetRasterBand(1)
areaBand = dataset.GetRasterBand(2)


nonzero = 0
total = 0
band = dataset.GetRasterBand(4)
noDataValue = band.GetNoDataValue()

array = band.ReadAsArray(0,0,band.XSize,band.YSize)
print numpy.where(array != noDataValue)
nonzero = numpy.where(array != noDataValue)[0].size
total = band.XSize*band.YSize

#for y in range(band.YSize):
#    scanline = band.ReadAsArray(0,y,band.XSize, 1)
#    nonzero += numpy.where(scanline != noDataValue)[0].size
#    total += band.XSize

#print 'RasterCount: ' + str(dataset.RasterCount)
#print 'upper left: (' + str(originX) + ', ' + str(originY) + ')'
#print 'pixel dimensions (' + str(pixelWidth) + ', ' + str(pixelHeight) + ')'
print 'nonzero elements: ' + str(nonzero) + ' total: ' + str(total) + ' percent: ' + str(float(nonzero)/total)
#print 'Time elapsed: ' + str(time.time()-startTime) + ' seconds'