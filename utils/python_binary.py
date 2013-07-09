import struct

FILE = './test.bin'
map = {}
map[(1,2)] = {0: [None, 0.25], 1: [.82, .45]}
map[(2,1)] = {0: [None, 0.25], 1: [.82, .45]}

def pickleBinaryMap(map,fileName):
    with open(fileName,'wb') as outFile:
        outFile.write(struct.pack('!i',len(map))) # ncoords
        for coord in map: 
            outFile.write(struct.pack('!h',coord[0])) # x coord
            outFile.write(struct.pack('!h',coord[1])) # y coord
            outFile.write(struct.pack('!b',len(map[coord]))) # n crops
            for cropId in map[coord]:
                outFile.write(struct.pack('!b',cropId)) # cropId
                for val in map[coord][cropId]: # crop value (either area or yield)
                    if val != None:
                        outFile.write(struct.pack('!f',val))
                    else:
                        outFile.write(struct.pack('!f',-1.0))

def unpickleBinaryMap(filename):
    with open(filename,'rb') as inFile:
        nCoords = struct.unpack('!i',inFile.read(4))[0] # ncoords
        #print nCoords
        globalMap = {}
        for i in range(nCoords):
            coord = (struct.unpack_from('!h',inFile.read(2))[0],struct.unpack_from('!h',inFile.read(2))[0]) # x, y coord
            #print coord
            nCrops = struct.unpack_from('!b',inFile.read(1))[0] #n crops
            #print nCrops
            cropMap = {}
            for j in range(nCrops):
                cropId = struct.unpack_from('!b',inFile.read(1))[0] #cropId
                #print cropId 
                cropArray = [None,None]
                for k in range(2):
                    cropVal = struct.unpack_from('!f',inFile.read(4))[0]  # crop value (either area or yield) 
                    if cropVal != -1.0:
                        cropArray[k] = cropVal
                cropMap[cropId] = cropArray
                #print cropArray
            globalMap[coord] = cropMap
        #print globalMap
        return globalMap

pickleBinaryMap(map,FILE)
print unpickleBinaryMap(FILE)
