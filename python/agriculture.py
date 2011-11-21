import os, sys, time, numpy
from osgeo import gdal
from osgeo import osr
from osgeo.gdalconst import *

FILENAME = '../sample_data/lulc_samp_cur'
CROPLANDID = 1

def die(message):
    print message
    sys.exit(1)

def getCoordsWithId(band,id):
    #Read dataset as 2D array, could kill memory -- should read as block
    array = band.ReadAsArray(0,0,band.XSize,band.YSize)
    return numpy.where(array == CROPLANDID)
    
    #Loop through dataset by hand, very slow
    #for x in range(0,band.XSize):
    #    for y in range(0,band.YSize):
    #        if array[y,x] == 1: 
    #            croplandCount += 1
    #print croplandCount

    #Loop through line by line, good balance between memory footprint and runtime
    #croplandCount = 0
    #for y in range(band.YSize):
    #    scanline = band.ReadAsArray(0,y,band.XSize, 1)
    #    croplandCount += numpy.where(scanline == CROPLANDID)[0].size
    #return croplandCount

#for timing runs
startTime = time.time()

#Register all drivers at once
gdal.AllRegister()

#open file
dataset = gdal.Open(FILENAME, GA_ReadOnly);

if dataset is None:
    die('Could not open ' + FILENAME)

if (dataset.RasterCount != 1):
    die('RasterCount should be 1 '+ dataset.RasterCount)

#Get georeference info, coordinates are for top left corners of pixels
geotransform = dataset.GetGeoTransform()
originX = geotransform[0]
originY = geotransform[3]
pixelWidth = geotransform[1]
pixelHeight = geotransform[5]

#define projection for indexing lat/lng
srs = osr.SpatialReference()
srs.ImportFromWkt(dataset.GetProjection())
srsLatLong = srs.CloneGeogCS() #gets the coordinate system
ct = osr.CoordinateTransformation(srs,srsLatLong) #transform from current projection to that projection's coordinate system


#Get first band, should be only one since it's LULC    
band = dataset.GetRasterBand(1)

#find the projected coordinates of all LULC with CROPLANDID
idCoords = getCoordsWithId(band,CROPLANDID)

#project to lat/lng
croplandLatLng = ct.TransformPoints(zip(idCoords[0],idCoords[1]))



#Results
print 'upper left: (' + str(originX) + ', ' + str(originY) + ')'
print 'pixel dimensions (' + str(pixelWidth) + ', ' + str(pixelHeight) + ')'
print 'X: ' + str(band.XSize) + ' Y: ' + str(band.YSize)
print 'Total cropland: ' + str(len(croplandLatLng)) + ' should be 21992'
print 'Time elapsed: ' + str(time.time()-startTime) + ' seconds'