from sets import Set
import os, sys
from numpy import *
try:
    from osgeo import ogr, gdal
    from osgeo.gdalconst import *
    import numpy
    use_numeric = False
except ImportError:
    import ogr, gdal
    from gdalconst import *
    import Numeric
class checkASCIIFile:
    def __init__(self, theFile):
        theFile = open(theFile, "r")
        self.line = theFile.readline()
        self.holder = []
        while self.line:
            self.holder.append(int(self.line.split(":")[0]))
            self.line = theFile.readline()
        theFile.close()

    def countDuplicatesInList(self, dupedList):
        uniqueSet = Set(item for item in dupedList)
        return [(item, dupedList.count(item)) for item in uniqueSet]
    
    def hasDuplicates(self):
        for x in self.countDuplicatesInList(self.holder): 
            if x[1] > 1:
                return True
        return False
    def isSorted(self):
        xhold = 0
        for x in self.holder:
            if xhold:
                if x < xhold:
                    return False
            xhold = x
        return True

def createfolders(theWorkspace):

    try:
        thefolders=["output","intermediate"]
        for folder in thefolders:
            if not os.path.exists(theWorkspace+"\\"+folder):
                os.mkdir(theWorkspace+"\\"+folder)
    except:
        print "Error creating folders"
        
def switchws(direction="up"):
    if direction=="down":
        gp.workspace=original_workspace+"\\intermediate"
    else:
        gp.workspace=original_workspace 


def convertRastertoArray(theDir, theRaster):
    use_numeric = True
    os.chdir(theDir) 
    # register all of the drivers
    gdal.AllRegister()
    
    # open the image
    ds = gdal.Open(theRaster, GA_ReadOnly)
    if ds is None:
        raise Exception, 'Could not open image'
 
    band = ds.GetRasterBand(1) # 1-based index
    # read data and add the value to the string
    thearray = band.ReadAsArray()
    return thearray 

def readFileIntoArray(theFile, skip=0):
    matrixfile = open(theFile)
    line = matrixfile.readline()
    if skip:
        line = matrixfile.readline()
    distmax = []
    while line:
        myarr = line.split()
        distmax.append(myarr)
        line = matrixfile.readline()
    matrixfile.close()
    return array(distmax).astype(float)

def convertArraytoRaster(theArray, outImage, theDir, thetrans, theproj):

    os.chdir(theDir)
    dims = theArray.shape 
    nx=dims[1] 
    ny=dims[0]
    nxarray=numpy.zeros(nx)  
    cal_factor=-83.0
    
    dst_format = 'HFA' 
    dst_datatype = gdal.GDT_Float32
    dst_options = []
    dst_file = outImage 
    dst_xsize = nx 
    dst_ysize = ny 
    dst_nbands = 1 
    
    
    geoTransform = thetrans
    proj = theproj
    
    driver = gdal.GetDriverByName( dst_format ) 
    dst_ds = driver.Create(dst_file, dst_xsize, dst_ysize, dst_nbands, 
    dst_datatype, dst_options) 
    
    
    for j in range(ny): 
        nxarray = theArray[j,:]
        nxarray.shape = (1,nx)
        dst_ds.GetRasterBand(1).WriteArray(nxarray,0,j) 
    
    dst_ds.SetGeoTransform(geoTransform)
    dst_ds.SetProjection(proj)
    
    dst_ds = None
    return

def getTransform (thedir, raster):
    os.chdir(thedir)
    ds = gdal.Open(raster, GA_ReadOnly)   
    geoTransform = ds.GetGeoTransform()
    ds = None
    return geoTransform
        
def getProj (thedir, raster):
    os.chdir(thedir)
    ds = gdal.Open(raster, GA_ReadOnly)   
    proj = ds.GetProjection()
    ds = None
    return proj

#function counts duplicates in a list
#arguments: python list
#return: list with count of dups
def countDuplicatesInList(dupedList):
    uniqueSet = Set(item for item in dupedList)
    return [(item, dupedList.count(item)) for item in uniqueSet]

#function checks if list has duplicates
#arguments: python list from countDuplicatesInList
#return: boolean, true is has duplicates
def ListHasDuplicates(theList):
    for x in countDuplicatesInList(theList): 
        if x[1] > 1:
            return True
    return False
#function sorts a text file before its used for reclassing raster
def sortreclassfile(thefile):
    newlines = []
    filetosort = open(thefile)
    lines = filetosort.readlines()
    for line in lines:
        splat = line.split(":")
        splat = map(int, splat)
        newlines.append(splat)
    filetosort.close()
    filetosort = open(thefile, "w")
    newlines = sorted(newlines, key=lambda cover: cover[0])
    for line in newlines:
        filetosort.write(str(line[0])+":"+str(line[1])+"\n")
    filetosort.close()
    
        
    





