# Import system modules
import string, arcpy, sys, os, shutil
from arcpy.sa import *
import numpy
import logging

#setup the logger
logging.basicConfig(format='%(asctime)s %(name)-18s %(levelname)-8s \
%(message)s', level=logging.DEBUG, datefmt='%m/%d/%Y %H:%M:%S ')
logger = logging.getLogger('scenario')

# Check out any necessary licenses
arcpy.CheckOutExtension("spatial")

# Set overwrite output
arcpy.env.overwriteOutput = 1

#global variables
outputDirectoryPrefix = ''
intermediateDirectoryPrefix = ''
probrasterDirectoryPrefix = ''
suitabilityDirectoryPrefix = ''
factorsDirectoryPrefix = ''
scenarioargs = {}
lulc_count = 0
factor_weight = 0
mainlcarray = [] #a multidimensional array with LULC class, COUNT, CHANGE (%) and No. pixels to change
onemask = numpy.array([]) #Create array used to mask pixels that have already been converted
landcoverarray = numpy.array([]) #array for holding landcover
masterarray = numpy.array([], dtype='S12')
masterarraysorted = numpy.array([], dtype='S12')
landcoverclassesarray = numpy.array([], dtype='int32')

debug = False #change to true to enter parameters manually



def execute(args):
    """Executes the scenario builder model that creates a new dataset based on scenario rules.
        scenarioargs - is an arguments dictionary with the following entries:
        scenarioargs['workspace_dir'] - the directory that will write output
            and other temporary files during calculation. (required)
        scenarioargs['scratch_dir'] - is scratch workspace
        scenarioargs['landcover'] - is the current land use/land cover map (required)
        scenarioargs['landattrib'] - is a DBF dataset with land attributes.
        returns nothing"""
    
    global probrasterDirectoryPrefix, outputDirectoryPrefix, probrasterDirectoryPrefix, suitabilityDirectoryPrefix, intermediateDirectoryPrefix, factorsDirectoryPrefix
    global scenarioargs
    
    if not debug:
        arcpy.env.workspace = arcpy.GetParameterAsText(0)
        lcraster = arcpy.GetParameterAsText(1)
        landattributes = arcpy.GetParameterAsText(2)
        overridelayer = arcpy.GetParameterAsText(3)
        constraintslayer = arcpy.GetParameterAsText(4)
        factorstable = arcpy.GetParameterAsText(5)
        factorweight = arcpy.GetParameterAsText(6)
        resolution = arcpy.GetParameterAsText(7)
        useprobability = arcpy.GetParameterAsText(8)
        usefactors = arcpy.GetParameterAsText(9)
        useproximity = arcpy.GetParameterAsText(10)
        resultOutput = arcpy.GetParameterAsText(11)
        donotrandomize = arcpy.GetParameterAsText(12)
    
        scenarioargs = {}
        scenarioargs['workspace_dir'] = arcpy.env.workspace + os.sep
        scenarioargs['landcover'] = str(lcraster)
        scenarioargs['landattrib'] = landattributes
        scenarioargs['override'] = overridelayer
        scenarioargs['constraints'] = constraintslayer
        scenarioargs['factors'] = factorstable
        scenarioargs['resolution'] = resolution
        if resolution:
            scenarioargs['resolution'] = int(str(resolution))
        else:
            scenarioargs['resolution'] = None
        scenarioargs['result'] = resultOutput
        scenarioargs['useprobability'] = True if useprobability == 'true' else False
        scenarioargs['usefactors'] = True if usefactors == 'true' else False
        scenarioargs['useproximity'] = True if useproximity == 'true' else False
        factor_weight = float(factorweight)
        scenarioargs['notrandomized'] = True if donotrandomize == 'true' else False
        if scenarioargs['notrandomized']:
            numpy.random.seed(42)
    else:
        
        #debugging
        scenarioargs['workspace_dir'] = 'C:\\DTL\\'
        scenarioargs['landcover'] = r'C:\DTL\Input\dtl.gdb\ludtl_rc'
        scenarioargs['landattrib'] = r'C:\DTL\Input\land_attributes.dbf'
        scenarioargs['override'] = ''#'C:\\scenarios_test\\input\\override.shp'
        scenarioargs['constraints'] = '' #'C:\\scenarios_test\\input\\constraints.shp'
        scenarioargs['factors'] = 'C:\\DTL\\Input\\factors.dbf'
        scenarioargs['resolution'] = ''
        if scenarioargs['resolution']:
            scenarioargs['resolution'] = int(str(scenarioargs['resolution']))
        else:
            scenarioargs['resolution'] = None
        scenarioargs['useprobability'] = True 
        scenarioargs['usefactors'] = True
        scenarioargs['useproximity'] = False   
        global factor_weight
        factor_weight = 0.5
        scenarioargs['result'] = ''
    
    overrideonly = False


    
    #ensure that either factors or probability or both are used
    if not scenarioargs['usefactors'] and not scenarioargs['useprobability'] and not scenarioargs['override']:
        raise Exception, 'You must use either converstion factors (eg roads) or probability of change or enter an override layer'
    
    #if only override given, do not convert pixels
    if scenarioargs['override'] and not scenarioargs['usefactors'] and not scenarioargs['useprobability']:
        overrideonly = True
    
    ##adjust resolution
    desc = arcpy.Describe(scenarioargs['landcover'])
    if scenarioargs['resolution'] > 0:
        if int(desc.meanCellHeight) <> scenarioargs['resolution']:
            
            #check if there is already a landcover with that resolution
            if arcpy.Exists(scenarioargs['workspace_dir']+'Intermediate'+os.sep+'lcresample'):
                desc = arcpy.Describe(scenarioargs['workspace_dir']+'Intermediate'+os.sep+'lcresample')
                if int(desc.meanCellHeight) <> scenarioargs['resolution']:
                    AddMsgAndPrint('Resampling landcover to ' + str(scenarioargs['resolution']))
                    arcpy.Resample_management(scenarioargs['landcover'], scenarioargs['workspace_dir']+'Intermediate'+os.sep+'lcresample', scenarioargs['resolution'])
            else:
                AddMsgAndPrint('Resampling landcover to '+scenarioargs['resolution'])
                arcpy.Resample_management(scenarioargs['landcover'], scenarioargs['workspace_dir']+'Intermediate'+os.sep+'lcresample', scenarioargs['resolution'])
            scenarioargs['landcover'] = str(scenarioargs['workspace_dir']+'Intermediate'+os.sep+'lcresample')
    
    #set arc environment
    arcpy.env.cellSize = scenarioargs['landcover']
    arcpy.env.extent = scenarioargs['landcover']
    

    
    prepare_landattrib_array()
    
    #These lines sets up the output directory structure for the workspace
    outputDirectoryPrefix = scenarioargs['workspace_dir'] + 'Output' + os.sep
    intermediateDirectoryPrefix = scenarioargs['workspace_dir'] + \
        'Intermediate' + os.sep
    probrasterDirectoryPrefix = scenarioargs['workspace_dir']  + 'Intermediate'+os.sep+'probability' + os.sep
    suitabilityDirectoryPrefix = scenarioargs['workspace_dir'] + 'Intermediate'+os.sep+'suitability' + os.sep
    factorsDirectoryPrefix = scenarioargs['workspace_dir'] + 'Intermediate'+os.sep+'factors' + os.sep
    scratchDirectoryPrefix = scenarioargs['workspace_dir'] + 'scratch' + os.sep
    
    for d in [outputDirectoryPrefix, intermediateDirectoryPrefix, suitabilityDirectoryPrefix, probrasterDirectoryPrefix, scratchDirectoryPrefix]:
        if not os.path.exists(d):
            AddMsgAndPrint('creating directory %s', d)
            os.makedirs(d)
    #set scratch directory
    arcpy.env.scratchWorkspace = scratchDirectoryPrefix

    #validate the input
    validationmsg = validate_input()
            
    if overrideonly:
        outRaster = apply_override(None, 1)
        #outSetNull = SetNull(intRas, intRas, "VALUE = 0")   
        outRaster.save(outputDirectoryPrefix+'result') 
        if scenarioargs['result']:
            outRaster.save(scenarioargs['result'])
            arcpy.DefineProjection_management(scenarioargs['result'], scenarioargs['landcover'])
        else:
            outRaster.save(outputDirectoryPrefix+'result')
            arcpy.DefineProjection_management(outputDirectoryPrefix+'result', scenarioargs['landcover'])        
        return    
    
    #create the landcover array
    global landcoverarray
    landcoverarray = arcpy.RasterToNumPyArray(scenarioargs['landcover'])
    


    #Create array used to mask pixels that have already been converted
    global onemask
    onemask = numpy.ones(landcoverarray.shape, dtype=numpy.int)

    #Mask out covers that are increasing so that they dont get converted    
    mask_increasing_areas()


    
    #prepare an array with the LULC numbers, Count of pixels, Change and number of pixels to change(goal)
    prepare_lulc_array()   
    
    #create probaility rasters
    if scenarioargs['useprobability']:
        create_probability_conversion_rasters()

    
    #combine factors
    if scenarioargs['usefactors']:
        combine_factors()
        
    #combine suitability
    combine_suitability_rasters()
        
    #Apply constraints to the suitability layers
    if scenarioargs['constraints']:
        apply_constraints()
    

    
    #Apply proximity to the suitability layers
    if scenarioargs['useproximity']:
        apply_proximity()
    
    #iterate though the cover types and change the pixels
    convert_pixels()
    
def mask_increasing_areas():
    """Creates a mask of covers that are growing so that they are not converted"""
    AddMsgAndPrint('Creating mask...')
    cur_increase = arcpy.SearchCursor(scenarioargs['landattrib'])
    reclasstable = open(scenarioargs['workspace_dir']+"increasing.asc","w")
    for lcclass in cur_increase:
        if lcclass.change < 0:
            reclasstable.write(str(int(lcclass.LULC))+": 0\n")
    reclasstable.close()
    outReclass = arcpy.sa.ReclassByASCIIFile(scenarioargs['landcover'],scenarioargs['workspace_dir']+"increasing.asc")
    outSetNull = SetNull(outReclass, outReclass, "VALUE > 0")
    outSetNullZeroOne = Con(IsNull(outSetNull), 0, 1)
    outSetNullZeroOne.save(scenarioargs['workspace_dir']+'Intermediate'+os.sep+"increasemask")
    
def create_probability_conversion_rasters():
    """Creates rasters of probability of conversion from the likelihood values in the land attributes table"""
    if not scenarioargs['useprobability']:
        return
    lulccodes = masterarray[:,0]
    lulcnames = masterarray[:,1]
    lulcchange = masterarray[:,3]
    lulcchange = lulcchange.astype(float)
    #suitabulity indices start from index 5 (0 based)
    AddMsgAndPrint('Creating probability rasters...')
    position = 0
    for x in lulccodes:
        reclassablecover = None
        if lulcchange[position] > 0:
            AddMsgAndPrint(lulcnames[position])
            reclasstable = open(intermediateDirectoryPrefix + 'reclass.asc', "w")
            values = masterarray[:,8 + position] #the probability values start at position 8(0 based) in the masterarray
            values = values.astype(float)
            maxvalue = max(values)
            if maxvalue > 10:
                adjusted = map(lambda p: int((p / maxvalue) * 10), values)
            else:
                adjusted = values
            for y in range(len(adjusted)):
                if not cover_type_exists(int(lulccodes[y])):
                    adjusted[y] = 0
                if adjusted[y] > 0:
                    reclassablecover = True
                reclasstable.write(lulccodes[y] + ':' + str(adjusted[y])+"\n")
            reclasstable.close()
            if reclassablecover:
                outReclass = ReclassByASCIIFile(scenarioargs['landcover'], intermediateDirectoryPrefix + 'reclass.asc')
                arcpy.Delete_management(intermediateDirectoryPrefix + 'reclass.asc')
                outReclass.save(intermediateDirectoryPrefix + 'temp')
                arcpy.env.extent = intermediateDirectoryPrefix + 'temp'
                outSetNull = SetNull(intermediateDirectoryPrefix + 'temp', intermediateDirectoryPrefix + 'temp', "VALUE = 0")
                outSetNull.save(probrasterDirectoryPrefix + lulcnames[position])
            #arcpy.Copy_management(intermediateDirectoryPrefix + 'probability' + os.sep + lulcnames[position], suitabilityDirectoryPrefix + lulcnames[position])
        position += 1
        
    
def apply_constraints():
    """Apply constraints to the suitability raster"""

    arcpy.PolygonToRaster_conversion(scenarioargs['constraints'], 'porosity', intermediateDirectoryPrefix+'constraints', '','',scenarioargs['landcover'])
    conRaster = Con(IsNull(Raster(intermediateDirectoryPrefix+'constraints')), 1, Raster(intermediateDirectoryPrefix+'constraints'))
    AddMsgAndPrint('Applying constraints...')
    for lc in mainlcarray:
        if lc[2] > 0:
            AddMsgAndPrint(lc[4])
            outTimes = Int(Times(conRaster, suitabilityDirectoryPrefix + (lc[4].lower())))
            outTimes.save(suitabilityDirectoryPrefix + (lc[4].lower()+'1'))
            arcpy.Delete_management(suitabilityDirectoryPrefix + (lc[4].lower()))
            arcpy.Copy_management(suitabilityDirectoryPrefix + (lc[4].lower()+'1'), suitabilityDirectoryPrefix + (lc[4].lower()))
            arcpy.Delete_management(suitabilityDirectoryPrefix + (lc[4].lower()+'1'))
            
def apply_proximity():
    """Apply proximity suitability to the suitability raster"""
    lulcchange = masterarray[:,3]
    lulcchange = lulcchange.astype(float) 
    arcpy.env.CellSize = scenarioargs['landcover']
    arcpy.env.extent = scenarioargs['landcover']   
    AddMsgAndPrint('Applying proximity...')
    for lc in masterarray:
        if float(lc[3]) > 0 and (int(lc[0]) in landcoverclassesarray) and lc[7] == 1:
            AddMsgAndPrint(lc[1])
            #create pixels with current cover, set everything else to null
            outSetNull = SetNull(scenarioargs['landcover'], 1, "VALUE <> "+lc[0])
            #create a distance raster from current cover
            outEucDistanceProx = EucDistance(outSetNull, float(lc[7]), arcpy.env.CellSize) 
            #find the max value in the distance raster
            arcpy.CalculateStatistics_management(outEucDistanceProx)
            maxvalprox = arcpy.GetRasterProperties_management(outEucDistanceProx, "MAXIMUM")
            InMaxDistance = int(float(maxvalprox.getOutput(0)))
            #normalize
            outNorm = (Int(10 - ((outEucDistanceProx / InMaxDistance) * 10))) / 10.0
            #change zero values to 1
            outProx = Con(outNorm, 1, outNorm, "VALUE = 0")
            #outProx.save(suitabilityDirectoryPrefix + (lc[1].lower()+'1'))
            nullProx = IsNull(outProx)  
            finalProx = Con(IsNull(outProx), 1, outProx)
            proximityraster = finalProx + Raster(suitabilityDirectoryPrefix + lc[1].lower())   
            #proximityraster.save(suitabilityDirectoryPrefix + lc[1].lower())

def get_pixel_count(num):
    """Get pixel count of given land cover class code (num)"""
    cur = arcpy.SearchCursor(scenarioargs['landcover'])    
    for lc in cur:
        if lc.VALUE == num:
            del cur
            return lc.COUNT
    del cur
    return 0

def cover_type_exists(lulccode):
    """Check if a given landcover type exists in landcover given land cover class code (num)"""
    cur = arcpy.SearchCursor(scenarioargs['landcover'])    
    for lc in cur:
        if lc.VALUE == lulccode:
            del cur
            return True
    del cur
    return False

def validate_input():
    """Check if the number of land cover types in the landcover raster
    matches the number of classes given in the land attribute file
    If there is a validation error, a message will be return otherwise nothing is returned
    """
    AddMsgAndPrint('Validating input...')
    validationresult = []
    cur = arcpy.SearchCursor(scenarioargs['landcover'])  
    temp = []
    for lc in cur:
        temp.append(lc.VALUE)
    lcrasterarray = numpy.array(temp) #array of land cover classes from the landcover raster
    global landcoverclassesarray
    landcoverclassesarray = lcrasterarray.astype('int32')
    #check that the landcover classes is unique
    if not check_array_unique(lcrasterarray):
        validationresult.append('Landcover raster classes are not unique')
        return validationresult
    cur = arcpy.SearchCursor(scenarioargs['landattrib']) 
    temp = []
    for lc in cur:
        temp.append(lc.LULC) 
    global lulc_count
    lulc_count = len(temp)
    lcattribarray = numpy.array(temp) #array of land cover classes from the landcover raster  

    #check that the land attrib classes is unique
    if not check_array_unique(lcrasterarray):
        validationresult.append('Land attribute classes are not unique')
        return validationresult
    
    cur = arcpy.SearchCursor(scenarioargs['landattrib']) 
    #check that the landcover classes and land attribute classes are the same
    if count_array_elements(numpy.intersect1d(lcrasterarray, lcattribarray, True)) != count_array_elements(lcattribarray):
        validationresult.append('Landcover raster and land attribute files do not match')
        return validationresult
    
    #check that the cover shortnames do not exceed 8 characters
    cur = arcpy.SearchCursor(scenarioargs['landattrib']) 
    for lc in cur:
        if len(lc.SHORTNME) > 8:
            raise Exception, 'The short names (Field SHORTNME) for your cover types should not be longer than 8 characters.  You entered '+ lc.SHORTNME +'. Check the table '+ scenarioargs['landattrib']

    cur = arcpy.SearchCursor(scenarioargs['landattrib']) 
    temp = []
    for lc in cur: 
        temp.append(lc.SHORTNME)
    if not check_array_unique(numpy.array(temp)):
        raise Exception, 'The short names (Field SHORTNME) for your cover types should be unique.  Check the table '+ scenarioargs['landattrib']
    
    return 0


    
def check_array_unique(arr):
    """Check if an array is unique. return false if not"""
    newarr = numpy.unique(arr)
    if count_array_elements(arr) == count_array_elements(newarr):
        return True
    return False
    
def count_array_elements(arraytocount):
    """Count number of elements in array"""
    count = 0
    for e in arraytocount:
        count += 1
    return count

def prepare_lulc_array():
    """Prepare an array with the LULC classes, Count, change values and number of 
    pixels to change for each class.  Assign this to a global variable"""
    global mainlcarray
    cur = arcpy.SearchCursor(scenarioargs['landattrib']) 
    out = []
    for lc in cur:
        lulcrow = []
        lulcrow.append(int(lc.LULC))
        count = get_pixel_count(lc.LULC)
        lulcrow.append(count)
        lulcrow.append(lc.CHANGE)
        lulcrow.append(int(count * lc.CHANGE * 0.01))
        lulcrow.append(lc.SHORTNME)
        out.append(lulcrow)
    mainlcarray = out
    return mainlcarray

def prepare_landattrib_array():
    """Prepare an array with all landcover data  the array has this structure
    0 LULC    The land cover code  
    1 SHORTNME  The land cover short name (max 12 characters) 
    2 COUNT  Number of pixels in each LULC class  
    3 %CHANGE  The percentage change expected on that cover  
    4 PIXELCHANGE  The number of pixels expected to change
    5 PRIORITY  The priority of this landcover (objective) 
    6 PROXIMITY Does this landcover suitability affected by proximity, eg do pixels
                closer to agriculture have higher chances of converting to agriculture [0,1] 
    7 PROXDIST  At what distance does the proximity influence die? (in meters, default 10,000m) 
    8 F1...Fn  The matching landcover probability score for the matrix 
    """
    global masterarray
    global masterarraysorted
    lccur = arcpy.SearchCursor(scenarioargs['landattrib'])
    
    #get the image size
    imageDesc = arcpy.Describe(scenarioargs['landcover'])
    
    lccodes = []
    for lc in lccur:
        lccodes.append(int(lc.LULC))
    lccur = arcpy.SearchCursor(scenarioargs['landattrib'])
    lulclist = []
    for lc in lccur:
        lulcrow = []
        lulcrow.append(int(lc.LULC))
        lulcrow.append(lc.SHORTNME)
        count = get_pixel_count(lc.LULC)
        lulcrow.append(count)
        lulcrow.append(lc.CHANGE)
        if count:
            #if the cover type already exists the change is based on the current extent
            lulcrow.append(int(count * lc.CHANGE * 0.01))
        else:
            #else the change is based on total size of landscape
            lulcrow.append(int(imageDesc.height * imageDesc.width * lc.CHANGE * 0.01))
        lulcrow.append(lc.PRIORITY)
        lulcrow.append(lc.PROXIMITY)
        lulcrow.append(lc.PROXDIST)
        for cov in lccodes:
            lulcrow.append(lc.getValue('F'+str(cov)))
        lulclist.append(lulcrow)
    lulclistsorted = sorted(lulclist, key=lambda cover: cover[5])
    masterarray = numpy.array(lulclist, dtype='S12')
    masterarraysorted = numpy.array(lulclistsorted, dtype='S12')
        
    

def random_select_elements_grow(a, numtoconvert, seedpc):
    """Randomly select pixels to convert and grow from the selected pixels"""
    x = a.shape[0]
    y = a.shape[1]
    #if not b.any():
    b = numpy.zeros(a.shape, dtype=int)
    c = numpy.where(a==1)
   
    numtoconvertseed = int(numtoconvert * seedpc/100.0)
    numtoconvertseed = 1 if numtoconvertseed < 1 else numtoconvertseed
    numtoconvertseed = 2 if numtoconvertseed == 1 else numtoconvertseed
    numtoconvertseed = 50 if numtoconvertseed > 50 else numtoconvertseed
    numtoconvertgrow = numtoconvert - numtoconvertseed
   
    candidates = numpy.where(a==1)
    pool = []
    #print len(c[0])
    for z in range(len(c[0])):
        coords = []
        coords.append(c[0][z])
        coords.append(c[1][z])
        pool.append(coords)
    #print pool
    for num in range(numtoconvertseed):
        #print len(pool)
        if len(pool):
            randomcoord = numpy.random.randint(0,len(pool))
            b[pool[randomcoord][0],pool[randomcoord][1]] = 1
            del pool[randomcoord]   
       
    for num in range(numtoconvertgrow):
        reset = 1
        #print len(pool)
        ceilfactor = 1000 #maximum number of randomizations it can do on each cell
        while reset and len(pool) and ceilfactor:
            randomcoord = numpy.random.randint(0,len(pool))
            xcoord = pool[randomcoord][0]
            ycoord = pool[randomcoord][1]
            xleft = xcoord if xcoord == 0 else xcoord-1
            xright = xcoord if xcoord >= (x-1) else xcoord+1
            ybottom = ycoord if ycoord == 0 else ycoord-1
            ytop = ycoord if ycoord >= (y-1) else ycoord+1
            if b[xleft][ycoord]==1 or b[xcoord][ytop]==1 or b[xright][ycoord]==1 or b[xcoord][ybottom]==1:
                b[xcoord,ycoord] = 1
                reset = 0
            ceilfactor -=1
        if len(pool):
            del pool[randomcoord]   
    return b
        
        
def convert_pixels():
    """Convert the pixels for each land cover type given its suitability raster"""
    global onemask
    global landcoverarray
    lulccodes = masterarraysorted[:,0]
    lulcnames = masterarraysorted[:,1]
    lulcchange = masterarraysorted[:,4]
    lulcchange = lulcchange.astype(float)
    position = 0
    AddMsgAndPrint('Converting pixels...')
    for x in lulccodes:
        if lulcchange[position] > 0:
            land_cover_shortname = lulcnames[position]
            number_pixels_tochange = int(lulcchange[position])
            land_cover_code = int(lulccodes[position])
            AddMsgAndPrint(land_cover_shortname)
            #if the cover is decreasing (change is -ve) skip to next cover
            if number_pixels_tochange <= 0:
                continue
            #read the suitability raster for this cover
            suitability_raster = suitabilityDirectoryPrefix + land_cover_shortname.lower()
            
            #apply increasemask to mask out increasing areas
            suitability_raster = Times(suitability_raster, Raster(intermediateDirectoryPrefix+'increasemask'))
    
            # Create a dictionary of count of pixels in suitability raster.
            
            suitability_raster_cur = arcpy.SearchCursor(suitability_raster)
            suit={}
            for suitclass in suitability_raster_cur:
                suit[suitclass.VALUE] = suitclass.COUNT    
            
            #create array of the suitability raster
            suitabilityarray = arcpy.RasterToNumPyArray(suitability_raster)
            suitabilityarray = numpy.around(suitabilityarray)
            suitabilityarray = suitabilityarray.astype('int32')
            newclass = int(land_cover_code)
            for x in range(10, 0, -1):
                if number_pixels_tochange > 0:
                    if x in suit:
                        if suit[x] < number_pixels_tochange:
                            #convert
                            criterion = (suitabilityarray==x) & onemask
                            #if number needed is less than number available, take random sample.
                           
        
                        else: #suitable pixels are more than the number required so pick at random
                            criterion = (suitabilityarray==x) & onemask
                            newarray = numpy.where(criterion > 1, 0, criterion)
                            if newarray.sum() < 1:
                                continue
                            ratio = number_pixels_tochange / float(suit[x])
                            if ratio > 0.5:
                                seed = 5
                            else:
                                seed = 1
                            if False:
                                cover.SHORTNME == r"Built"
                                orig = arcpy.RasterToNumPyArray("orig")
                                origval = numpy.where(orig>1, 0, orig)
                                criterion = random_select_elements_grow(criterion, number_pixels_tochange, seed, origval)
                            else:
                                criterion = random_select_elements_grow(criterion, number_pixels_tochange, seed)
                        #change pixel values
                        landcoverarray = numpy.where(criterion, newclass, landcoverarray)
                        onemask = numpy.where(criterion, 0, onemask)
                        number_pixels_tochange -= suit[x]
                        if number_pixels_tochange < 0: number_pixels_tochange = 0  
                else:
                    break;
        position += 1
    # Convert the landcover array back to raster
    x = arcpy.GetRasterProperties_management(scenarioargs['landcover'], "LEFT")
    y = arcpy.GetRasterProperties_management(scenarioargs['landcover'], "BOTTOM")
    origin = arcpy.Point(x.getOutput(0),y.getOutput(0))
    #apply the override
    if scenarioargs['override']:
        landcoverarray = apply_override(landcoverarray)
    
    newRaster = arcpy.NumPyArrayToRaster(landcoverarray, origin, scenarioargs['landcover'])
    intRas = Int(newRaster)
    outSetNull = SetNull(intRas, intRas, "VALUE = 0")   
    if scenarioargs['result']:
        outSetNull.save(scenarioargs['result'])
        arcpy.DefineProjection_management(scenarioargs['result'], scenarioargs['landcover'])
    else:
        outSetNull.save(outputDirectoryPrefix+'result')
        arcpy.DefineProjection_management(outputDirectoryPrefix+'result', scenarioargs['landcover'])        

    
def apply_override(arr, overrideonly=False):
    """Explicitly change the landcover based on given polygons(override)"""
    arcpy.env.mask = scenarioargs['landcover']
    if scenarioargs['override']:
        arcpy.PolygonToRaster_conversion(scenarioargs['override'], 'newclass', intermediateDirectoryPrefix+'override', '','',scenarioargs['landcover'])
        conRaster = Con(IsNull(Raster(intermediateDirectoryPrefix+'override')), 0, Raster(intermediateDirectoryPrefix+'override'))
        if overrideonly:
            return Con(conRaster, conRaster, Raster(scenarioargs['landcover']))
        else:
            overridearray = arcpy.RasterToNumPyArray(conRaster)
            output = numpy.where(overridearray, overridearray, arr)
            return output
    return

def combine_factors():
    if not scenarioargs['usefactors']:
        return 
    if not arcpy.Exists(scenarioargs['factors']):
        return
    

    #check that the cover shortnames in factors table do not exceed 8 characters
    cur = arcpy.SearchCursor(scenarioargs['factors']) 
    #for lc in cur:
        #if len(lc.COVER) > 8:
            #pass            
            #raise Exception, 'The cover names (Field COVER) for your cover types should not be longer than 8 characters.  You entered '+ lc.COVER +'. Check the table '+ scenarioargs['factors']
     
    lulccodes = masterarray[:,0]
    lulcnames = masterarray[:,1]
    lulcchange = masterarray[:,3]
    lulcchange = lulcchange.astype(float)
    #suitabulity indices start from index 5 (0 based)
    position = 0
    #step through each cover from the land attributes table 
   # AddMsgAndPrint('Creating factor rasters...')
    for x in lulccodes:
        if lulcchange[position] <= 0:
            position += 1
            continue
        
        factors = scenarioargs['factors']

        #Process the factors for this land cover class.
        #If factors table given for this cover, process them else skip
        if arcpy.Exists(factors):
            AddMsgAndPrint(lulcnames[position])
            cur_factors = arcpy.SearchCursor(factors)
            # Check if the data for the factor exists.
            for factor in cur_factors:
                if not arcpy.Exists(factorsDirectoryPrefix+factor.getValue("cover") +"_"+ factor.getValue("factorname")[:3]):
                    raise Exception, "Dataset "+ factor.getValue("cover")+"_"+factor.getValue("factorname")+" does not exist."
            cur_factors = arcpy.SearchCursor(factors)
            coverfactors = []
            count = 0
            out = None
            for changefactor in cur_factors:
                if changefactor.cover.lower() == lulcnames[position].lower():
                    AddMsgAndPrint(changefactor.factorname)
                    fact = Con(IsNull(factorsDirectoryPrefix + changefactor.cover + "_" + changefactor.factorname[:3]),0,factorsDirectoryPrefix + changefactor.cover + "_" + changefactor.factorname[:3])
                    if count > 0:
                        out = out + Times(fact, (changefactor.weight))
                    else:
                        out = Times(fact, (changefactor.weight))
                    count += 1    
            if arcpy.Exists(out) or out:
                out.save(factorsDirectoryPrefix + lulcnames[position])
        position += 1
        
def combine_suitability_rasters():
    AddMsgAndPrint('Combining suitability rasters...')
    lulccodes = masterarray[:,0]
    lulcnames = masterarray[:,1]
    lulcchange = masterarray[:,3]
    lulcchange = lulcchange.astype(float)
    #suitabulity indices start from index 5 (0 based)
    position = 0
    #step through each cover from the land attributes table 
    
    #if both factorfile and probablity matrix given, combine
    AddMsgAndPrint(scenarioargs['usefactors'])
    AddMsgAndPrint(scenarioargs['useprobability'])

    if scenarioargs['usefactors'] and scenarioargs['useprobability']:
        for x in lulccodes:
            if lulcchange[position] <= 0:
                position += 1
                continue
            
            if arcpy.Exists(probrasterDirectoryPrefix + lulcnames[position]) and arcpy.Exists(factorsDirectoryPrefix + lulcnames[position]):
                suitraster = Int(Times((probrasterDirectoryPrefix + lulcnames[position]), (1 - factor_weight)) + Times((factorsDirectoryPrefix + lulcnames[position]), factor_weight))
            elif arcpy.Exists(probrasterDirectoryPrefix + lulcnames[position]):
                suitraster = Int(Raster(probrasterDirectoryPrefix + lulcnames[position]))
            elif arcpy.Exists(factorsDirectoryPrefix + lulcnames[position]):
                suitraster = Int(Raster(factorsDirectoryPrefix + lulcnames[position]))
            suitraster.save(suitabilityDirectoryPrefix + lulcnames[position])
            position += 1
    elif scenarioargs['usefactors']:
        for x in lulccodes:
            if lulcchange[position] <= 0:
                position += 1
                continue
            if arcpy.Exists(factorsDirectoryPrefix + lulcnames[position]):
                suitraster = Int(Raster(factorsDirectoryPrefix + lulcnames[position]))
            suitraster.save(suitabilityDirectoryPrefix + lulcnames[position])
            position += 1 
    elif scenarioargs['useprobability']:
        for x in lulccodes:
            if lulcchange[position] <= 0:
                position += 1
                continue
            elif arcpy.Exists(probrasterDirectoryPrefix + lulcnames[position]):
                suitraster = Int(Raster(probrasterDirectoryPrefix + lulcnames[position]))
            suitraster.save(suitabilityDirectoryPrefix + lulcnames[position])
            position += 1          
        
def AddMsgAndPrint(msg, severity=0):
    # Adds a Message (in case this is run as a tool)
    # and also prints the message to the screen (standard output)
    # 
    logger.debug(msg)

    # Split the message on \n first, so that if it's multiple lines, 
    #  a GPMessage will be added for each line
    try:
        for string in msg.split('\n'):
            # Add appropriate geoprocessing message 
            #
            if severity == 0:
                arcpy.AddMessage(string)
            elif severity == 1:
                arcpy.AddWarning(string)
            elif severity == 2:
                arcpy.AddError(string)
    except:
        pass
    
if __name__ == '__main__':
    args = sys.argv;
    execute(args)
    logger.info('Complete')