# Import system modules
import string, arcpy, sys, os, dbfpy
from arcpy.sa import *
import numpy
import logging

#setup the logger
logging.basicConfig(format='%(asctime)s %(name)-18s %(levelname)-8s \
%(message)s', level=logging.DEBUG, datefmt='%m/%d/%Y %H:%M:%S ')
logger = logging.getLogger('scenario_factors')

# Check out any necessary licenses
arcpy.CheckOutExtension("spatial")

# Set overwrite output
arcpy.env.overwriteOutput = 1

scenarioargs = {}
masterarray = numpy.array([])

scenarioargs['workspace_dir'] = arcpy.GetParameterAsText(0)#'C:\\DTL'

scenarioargs['landcover'] = arcpy.GetParameterAsText(1)#'C:\\DTL\\Input\\dtl.gdb\\ludtl_rc'
scenarioargs['factorfiles'] = arcpy.GetParameterAsText(2)#'C:\\DTL\\Input\\factordtl.dbf'
outputDirectoryPrefix = scenarioargs['workspace_dir']+ os.sep + 'Output' + os.sep
inputDirectoryPrefix = scenarioargs['workspace_dir'] + os.sep + \
    'Input' + os.sep
intermediateDirectoryPrefix = scenarioargs['workspace_dir'] + os.sep + \
    'Intermediate' + os.sep
probrasterDirectoryPrefix = scenarioargs['workspace_dir'] + os.sep + 'Intermediate'+os.sep+'probability' + os.sep
suitabilityDirectoryPrefix = scenarioargs['workspace_dir'] + os.sep + 'Intermediate'+os.sep+'suitability' + os.sep
factorsDirectoryPrefix = scenarioargs['workspace_dir'] + os.sep + 'Intermediate'+os.sep+'factors' + os.sep
scratchDirectoryPrefix = scenarioargs['workspace_dir'] + os.sep + 'scratch' + os.sep

arcpy.env.cellSize = scenarioargs['landcover']
arcpy.env.extent = scenarioargs['landcover'] 

#check that the cover shortnames do not exceed 8 characters
cur = arcpy.SearchCursor(scenarioargs['factorfiles']) 
for lc in cur:
    if len(lc.COVER) > 8:
        raise Exception, 'The cover names (Field COVER) should match the shortname used in the land attribute table and not be longer than 8 characters.  You entered '+ lc.COVER +'. Check the table '+ scenarioargs['factorfiles']

def execute(args):
    create_folders()
    arcpy.env.workspace = 'C:\\DTL'
    arcpy.env.scratchWorkspace = scenarioargs['workspace_dir']+os.sep+"scratch"
    createfactors()
def checktype(dataset):
    desc = arcpy.Describe(dataset)
    #print desc.datatype
    if desc.datatype == 'ShapeFile' or desc.datatype == 'FeatureClass':
        return desc.shapetype
    return desc.datatype
        
def createfactors():
    AddMsgAndPrint('Creating factor files')
    cur = arcpy.SearchCursor(scenarioargs['factorfiles'])
    for factor in cur:
        print inputDirectoryPrefix+factor.layer
        if arcpy.Exists(inputDirectoryPrefix+factor.layer):
            layertype = checktype(inputDirectoryPrefix+factor.layer)
            factorname = factor.factorname
            suitability = factor.suitabilit
            InMaxDistance = factor.dist
            AddMsgAndPrint(factorname + " on " + factor.cover)
            InSourceData = inputDirectoryPrefix+factor.layer
            OutDistanceRaster = intermediateDirectoryPrefix+"distgrid"
            if layertype=='Point' or layertype=='Polyline':

                
                #IncellSize = "30"
                IncellSize = arcpy.env.cellSize
                # Create distance raster showing suitability varying from the factor.
                
                outEucDistance = EucDistance(InSourceData, InMaxDistance, arcpy.env.cellSize)
                outEucDistance.save(OutDistanceRaster)
               
                # Flip the values based on the suitability values
                if suitability.lower() == "near":
                    val = 1
                    #gp.SingleOutputMapAlgebra_sa("10 - (int(distgrid "+"/ "+str(InMaxDistance)+" * 10))", "polytemp")
                    outRas = Int(10 - ((Raster(OutDistanceRaster) / InMaxDistance) * 10))
                    outRas.save(intermediateDirectoryPrefix+"polytemp")
                elif suitability.lower() == "far":
                    val = 10
                    #gp.SingleOutputMapAlgebra_sa("int(distgrid "+"/ "+str(InMaxDistance)+" * 10)", "polytemp")
                    outRas = (Int(Raster(intermediateDirectoryPrefix+"distgrid") / InMaxDistance) * 10)
                    outRas.save(intermediateDirectoryPrefix+"polytemp")
                else:
                    raise Exception, "Suitability score missing"
                # Change zero values to 1 for minimal chance.
                #gp.SingleOutputMapAlgebra_sa("Con(polytemp == 0, 1, polytemp)", "polytemp1")
                OutRas = Con(intermediateDirectoryPrefix+"polytemp", 1, intermediateDirectoryPrefix+"polytemp", "VALUE = 0")
                outRas.save(intermediateDirectoryPrefix+"polytemp1")
                # Fill the NoData areas.
                #gp.SingleOutputMapAlgebra_sa("Con(IsNull(polytemp1), "+val+", polytemp1)", factorname)
                
                OutNullRas = IsNull(intermediateDirectoryPrefix+"polytemp1")
                OutRas = Con(OutNullRas, val, intermediateDirectoryPrefix+"polytemp1", "VALUE = 1")
                outSetNull = SetNull(outRas == 0, outRas)
                outSetNull.save(factorsDirectoryPrefix+factor.cover.lower()+"_"+factorname[:3])
               
    
            elif layertype=='Polygon':
                dist = InMaxDistance
                OutRaster = factorsDirectoryPrefix+factor.cover.lower()+"_"+factorname[:3]
                #create a raster using the suitability field
                
                # Process: FeatureToRaster_conversion
                arcpy.FeatureToRaster_conversion(InSourceData, factor.suitfield, OutRaster, arcpy.env.cellSize)                

                #result.save(factorsfolder+factor.cover.lower()+"_"+factorname[:3])

    return

def create_folders():
    #These lines sets up the output directory structure for the workspace
    for d in [intermediateDirectoryPrefix, factorsDirectoryPrefix, scratchDirectoryPrefix]:
        if not os.path.exists(d):
            AddMsgAndPrint('creating directory %s', d)
            os.makedirs(d)
            
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
    AddMsgAndPrint('Complete')