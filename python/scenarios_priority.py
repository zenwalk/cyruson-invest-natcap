import scenarios_ahp as ahp
import numpy, arcpy
from dbfpy import dbf
import logging, os

#setup the logger
logging.basicConfig(format='%(asctime)s %(name)-18s %(levelname)-8s \
%(message)s', level=logging.DEBUG, datefmt='%m/%d/%Y %H:%M:%S ')
logger = logging.getLogger('scenario')

def execute():
    prioritymatrix = arcpy.GetParameterAsText(0) #'C:/scenarios_test/Input/ahptest.dbf' 
    landattrib = arcpy.GetParameterAsText(1) #'C:/scenarios_test/Input/land_attributes.dbf'
    if os.path.exists(prioritymatrix):
        dbpriority = dbf.Dbf(prioritymatrix)
    else:
        raise Exception, 'Priority matrix table required'
    records = dbpriority.recordCount
    fields = len(dbpriority.fieldDefs)
    priorityids = []
    for rec in dbpriority:
        priorityids.append(rec[0])
    
    #check that the land attributes table exists and create the db object
    if os.path.exists(landattrib):
        landattribdb = dbf.Dbf(landattrib)
        recordslandattrib = landattribdb.recordCount
        
        if records != recordslandattrib:
            raise Exception, 'The priority matrix and the land attributes matrix must have the same number of records'
        
        landattribids = []
        for rec in landattribdb:
            landattribids.append(rec[0])
        if priorityids != landattribids:
            raise Exception, 'The priority matrix and the land attributes matrix must have the same identifiers in the first column'
            
    
    #ensure that the matrix is of the correct shape
    if abs(records - fields) != 3:
        raise Exception, 'The matrix is not of the required shape.  It should only have 2 fields more than the number of records'
    
    if 'PRIORITY' in dbpriority.fieldNames:
        priorityfield = 'PRIORITY'
    elif 'priority' in dbpriority.fieldNames:
        priorityfield = 'priority'
    elif 'Priority' in dbpriority.fieldNames:
        priorityfield = 'Priority'
    else:
        raise Exception, 'PRIORITY field required in priority table'
    vals = []
    for rec in dbpriority:
        valrow = []
        for x in range(2,records+2):
            if not rec[x]:
                valrow.append(0.0)
            else:
                valrow.append(float(rec[x]))
        vals.append(valrow)
    
    arr = numpy.array(vals)
    
    for row in range(records):
        for col in range(records):
            if arr[row,col] == 0 or not arr[row,col]:
                arr[row,col] = 1./arr[col,row]

    priorityvalues = ahp.calculateWeights(arr,4)
    prioritydict = {}
    recno = 0
    for rec in dbpriority:
        rec[priorityfield] = priorityvalues[recno]
        if landattrib:
            prioritydict[int(rec['LULC'])] = priorityvalues[recno]
        recno += 1
        rec.store()
    dbpriority.close()   

    if landattrib:
        landattribdb = dbf.Dbf(landattrib)
        for rec in landattribdb:
            rec['priority'] = prioritydict[int(rec['LULC'])]
            rec.store()
        landattribdb.close()

    

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