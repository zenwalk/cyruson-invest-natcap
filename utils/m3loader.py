#The purpose of this script is to unpickle the sparse M3 crop data and convert
#into a single coherent data structure that allows indexes by coordinate
#to map to a map of cropids which in turn map to average yields and fractional
#land usage.

import os, sys, time, numpy, glob, re, struct, mmap

class M3Loader:
    # This dictionary maps coordinates to file positions, or crop data if
    # those file positions have already been looked up
    #_globalMap = {}
    #_inFileMap = None

    def __init__(self, FILENAME):
        self._globalMap = {}
        with open(FILENAME, 'r+b') as inFile:
            self._inFileMap = mmap.mmap(inFile.fileno(), 0)
            self._unpickleCoordinates(self._inFileMap)
            #print len(self._globalMap)
    
    def _unpickleCoordinates(self, file):
        """Looks up coordinates and file positions for later access, this way
        we don't load the entire dataset at once, only as needed"""
        
        def read(format, file):
             """Simplifies reading struct formatted strings from file"""
             return struct.unpack(format, file.read(struct.calcsize(format)))[0]
        #load cropid map
        self._cropIds = {}
        nCropIds = read('!B',file) # number of crop ids
        for i in range(nCropIds):
            cropNameLen = read('!B',file) #length of the crop name 
            cropName = read('!'+str(cropNameLen)+'s',file) #the crop name
            cropId = read('!B',file) #the id associated with that name
            self._cropIds[cropName] = cropId
            print cropName, cropId
        
        #load geotransform
        #read lngOrig, latOrig, lngWidth, latWidth  - from http://www.gdal.org/gdal_tutorial.html
        self._geoTransform = {}
        for v in [0, 3, 1, 5]:
            self._geoTransform[v] = read('!d', file)
            print self._geoTransform[v]
        
        nCoords = read('!i',file) # ncoords
        self._globalMap = {}
        latMax = 0
        lngMax = 0
        for i in range(nCoords):
            if i % 100000 == 0: #follow progress
                percent = i / float(nCoords)
                print percent
                sys.stdout.flush()
            coord = (read('!h', file), read('!h', file)) # latIndex, lngIndex
            if coord[0] > latMax:
                latMax = coord[0]
            if coord[1] > lngMax:
                lngMax = coord[1]
            nCrops = read('!B', file) #n crops
            #the -1 comes from the fact that we already read the nCrops, but we want to re-read it on the load
            self._globalMap[coord] = file.tell() - 1
            #jump over the crop data 
            file.seek((1 + 4 + 4) * nCrops, os.SEEK_CUR)
    
    def _getCropData(self, filePos):
        """Parses out the crop data for the given file position"""
        self._inFileMap.seek(filePos)
        nCrops = struct.unpack_from('!B', self._inFileMap.read(1))[0] #n crops
        
        cropMap = {}
        for j in range(nCrops):
            cropId = struct.unpack_from('!B', self._inFileMap.read(1))[0] #cropId
            cropArray = [None, None]
            for k in range(2):
                cropVal = struct.unpack_from('!f', self._inFileMap.read(4))[0]  # crop value (either area or yield) 
                if cropVal != -1.0:
                    cropArray[k] = cropVal
            cropMap[cropId] = cropArray
        return cropMap
    
    def getCropData(self, coord):
        """Report crop map index, if not cached, look it up
        
        coord -- a (latIndex,lngIndex) index into the global crop map"""
        if coord not in self._globalMap:
            return None
        if type(self._globalMap[coord]) is int:
            self._globalMap[coord] = self._getCropData(self._globalMap[coord])
        return self._globalMap[coord]
        
        
    def latLngToCoord(self, lat, lng):
        originLng = self._geoTransform[0]
        originLat = self._geoTransform[3]
        dimLng = self._geoTransform[1]
        dimLat = self._geoTransform[5]
    
        latIndex = int((lat - originLat) / dimLat)
        lngIndex = int((lng - originLng) / dimLng)
        return (latIndex, lngIndex)
    
if __name__ == '__main__':
    m3loader = M3Loader('../tmp_data/m3pickle.bin')        

    #in Colombia, should have some banana data
    lat = 8.083333
    lng = -67.74999

    cropData = m3loader.getCropData(m3loader.latLngToCoord(lat, lng))[m3loader._cropIds['banana']]
    print cropData, ' fractional area should equal ', 0.000780
    
        
    #globalMap is loaded, for fun look up some crop data    
#    for coord in globalMap:
#        print coord, globalMap[coord], getCropData(inFileMap, globalMap[coord])
#        cropData = getCropData(inFileMap, globalMap[coord])
#        for key in cropData:        
#            print cropData[key][0]
#            sys.stdout.flush()
#        sys.exit(0)

