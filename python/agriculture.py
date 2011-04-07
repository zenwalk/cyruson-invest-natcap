import numpy
import gdal
from gdalconst import *

#open file
FILENAME = '../sample_data/landuse_90'
dataset = gdal.Open(FILENAME, GA_ReadOnly);
CROPLANDID = 1

if dataset is None:
    print 'barf'

if (dataset.RasterCount != 1):
    print 'RasterCount should be 1'+ dataset.RasterCount
    
band = dataset.GetRasterBand(1) 


croplandCount = 0

print 'X: ' + str(band.XSize) + ' Y: ' + str(band.YSize) 

for i in range(band.YSize):
    scanline = band.ReadAsArray(0,i,band.XSize, 1, band.XSize, 1)
    croplandCount += numpy.where(scanline == CROPLANDID)[0].size

#calculate number of nan cells?        
print 'Total cropland: ' + str(croplandCount) + ' should be 9940'
