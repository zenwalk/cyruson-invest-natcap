import scenarios_ahp as ahp
import numpy, arcpy
from dbfpy import dbf
import logging, os

#setup the logger
logging.basicConfig(format='%(asctime)s %(name)-18s %(levelname)-8s \
%(message)s', level=logging.DEBUG, datefmt='%m/%d/%Y %H:%M:%S ')
logger = logging.getLogger('scenario')

def execute():
    scenarioargs = {}
    scenarioargs['workspace_dir'] = arcpy.GetParameterAsText(0)# 'C:\\scenarios_test'
    scenarioargs['landcover'] = arcpy.GetParameterAsText(1)#'C:\\scenarios_test\\Input\\landcover'
    landattrib = arcpy.GetParameterAsText(2)#'C:/scenarios_test/Input/land_attributes.dbf'
    
    lccur = arcpy.SearchCursor(scenarioargs['landcover'])
    
    desc = arcpy.Describe(scenarioargs['landcover'])
    cellsize = desc.meanCellHeight
    cellarea = cellsize / 10000.0
    
    pixels = {}
    for lc in lccur:
        pixels[lc.VALUE] = lc.COUNT
    del lccur
    

        
    
    landattribcur = arcpy.SearchCursor(landattrib)
    AddMsgAndPrint('LULC\tNAME\t\t\t% CHANGE\t\tCURRENT AREA(HA)\t\tNEW AREA(HA)')
    sumarea = 0
    sumnewarea = 0
    for rec in landattribcur:
        sumarea += (getcount(rec.LULC, pixels) * cellarea) 
        sumnewarea += (getcount(rec.LULC, pixels) * cellarea * (1 + (rec.CHANGE/100.0)))
        AddMsgAndPrint(str(int(rec.LULC))+'\t'+rec.NAME+'\t\t'+str(rec.CHANGE)+'\t\t'+str(getcount(rec.LULC, pixels) * cellarea)+'\t\t'+str(getcount(rec.LULC, pixels) * cellarea * (1 + (rec.CHANGE/100.0))))
    AddMsgAndPrint('Current total area:\t'+str(sumarea)+' Ha')
    AddMsgAndPrint('New total area:\t'+str(sumnewarea) +' Ha')
    AddMsgAndPrint('Net:\t'+str(sumarea - sumnewarea) +' Ha ('+ str(round(abs((sumarea - sumnewarea) / sumarea * 100),0))+'%)' )


def getcount(val, pixdict):
    if val in pixdict.keys():
        return pixdict[val]
    else:
        return 0
    

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
    execute()
    logger.info('Complete')