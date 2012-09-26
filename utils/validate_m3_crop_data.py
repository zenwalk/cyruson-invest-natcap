#The purpose of this script is to explore various values of data in the M3 crop 
#dataset

import os, sys, time, numpy, glob, re
from osgeo import gdal
from osgeo import osr
from osgeo.gdalconst import *

#Global constants

MAXBANDS = 2
PATH = '../tmp_data/175crops/'

#Functions

def die(message):
    """Error handler"""
    print message
    sys.exit(1)

def validateM3Data(band, cropName, outOfRange):
    """Test if any data seems out of range"""
    noDataValue = band.GetNoDataValue()
    array2d = band.ReadAsArray(0, 0, band.XSize, band.YSize)
    nonzeroValidIndices = numpy.where((array2d != noDataValue) & (array2d != 0))
    values = array2d[nonzeroValidIndices]
    for i in range(len(values)):
        coord = (nonzeroValidIndices[0][i], nonzeroValidIndices[1][i])
        if values[i] > outOfRange:
            print values[i], cropName, coord

#Register all drivers at once
gdal.AllRegister()

#regular expression to break the cropname out of the path 
#example (apple) from /home/joe/path/apple_5min.nc
cropNameRe = re.compile("^(.*/)([^/]*)_5min")
cropIds = {}

filenames = glob.glob(os.path.join(PATH, '*.nc'))

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
    
    harvestedAreaBand = dataset.GetRasterBand(1)
    yieldBand = dataset.GetRasterBand(2)
    
    #validateM3Data(yieldBand, cropName, 1e1)
    validateM3Data(harvestedAreaBand, cropName, 0.01)
    
    currentTime = time.time()-startTime
    totalTime += currentTime
    
    print ' time: '+ str(currentTime) 
    sys.stdout.flush()
print "totalTime , " + str(totalTime)
sys.stdout.flush()
